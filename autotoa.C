/* autotoa.C
 *
 * Iteratively determine template and TOAs from data.
 * Paul Demorest, 2008.
 */
#include "Pulsar/psrchive.h"
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/IntegrationExpert.h"
#include "Pulsar/BasicIntegration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/PolnProfile.h"
#include "Pulsar/ProfileShiftFit.h"
#include "Pulsar/PhaseWeight.h"
#include "Pulsar/WaveletSmooth.h"
#include "Pulsar/AdaptiveSmooth.h"
#include "MEAL/ScaledVonMises.h"
#include "model_profile.h"
#include "MJD.h"

#include "Error.h"
#include "dirutil.h"
#include "strutil.h"

#include <fstream>
#include <iostream>

#include <math.h>
#include <string.h>
#include <unistd.h>

using namespace std;
using namespace Pulsar;

void usage() {
    cout << "autotoa - Iteratively determine template and TOAs\n"
        "Usage: autotoa [options] file1 (file2 ..)\n"
        "Options:\n"
        "  -h       Help message\n"
        "  -v       Verbose output\n"
        "  -M       File names in metafile\n"
        "  -F       Fscrunch before timing\n"
        "  -I       Use invariant interval\n"
        "  -i num   Maximum iterations\n"
        "  -n num   Number of harmonics to use\n"
        "  -s file  Use initial standard/template given in file\n"
        "  -g width Use a single Gaussian of given width initially\n"
        "  -t file  Output TOAs to file\n"
        "  -S file  Output final template to file\n"
        << endl;
}

int main (int argc, char *argv[]) try {

    /* Process cmd line */
    int opt=0;
    bool verbose=false;
    bool use_tmpl=false;
    bool fscrunch=false;
    bool tscrunch=false;
    bool use_invint=false;
    char *metafile=NULL;
    int max_it=25, nharm=256;
    float snr_cutoff=1.00;
    float gauss_width=0.0;
    string tmplfile, toaoutfile="", tmploutfile="";
    while ((opt=getopt(argc,argv,"hvs:FTM:i:n:It:S:g:"))!=-1) {
        switch (opt) {
            case 's':
                use_tmpl=true;
                tmplfile = optarg;
                break;
            case 'F':
                fscrunch=true;
                break;
            case 'T':
                tscrunch=true;
                break;
            case 'M':
                metafile=optarg;
                break;
            case 'i':
                max_it = atoi(optarg);
                break;
            case 'v':
                Archive::set_verbosity(2);
                verbose=true;
                break;
            case 'n':
                nharm = atoi(optarg);
                break;
            case 'I':
                use_invint=true;
                break;
            case 't':
                toaoutfile = optarg;
                break;
            case 'S':
                tmploutfile = optarg;
                break;
            case 'g':
                gauss_width = atof(optarg);
                break;
            case 'h':
            default:
                usage();
                return(0);
                break;
        }
    }

    /* Warn if no outputs */
    if (toaoutfile=="" && tmploutfile=="")
        cerr << "Warning: No output will be produced! Use -t and/or -S." 
            << endl;

    /* Get list of files */
    vector<string> archives;
    if (metafile)
        stringfload(&archives, metafile);
    else
        for (int ai=optind; ai<argc; ai++) dirglob(&archives, argv[ai]);
    if (archives.empty()) {
        cerr << "No files given." << endl;
        usage();
        return(-1);
    }

    /* Fancy pointers */
    Reference::To<Archive> arch;
    Reference::To<Archive> tmplarch;
    Reference::To<Integration> subint;
    Reference::To<Profile> prof;
    Reference::To<Profile> tmplprof;
    std::vector< Reference::To<Profile> > sum;
    sum.resize(4);
    Reference::To<Profile> isum;
    Reference::To<PolnProfile> pprof;

    // Template smoothing
    //WaveletSmooth smooth;
    AdaptiveSmooth smooth;

    /* Load initial template */
    if (use_tmpl) try {
        tmplarch = Archive::load(tmplfile);
        /* TODO: allow freq-dependent template? */
        tmplarch->fscrunch();
        tmplarch->tscrunch();
        if (use_invint)
            //tmplarch->convert_state(Signal::Invariant);
            tmplarch->convert_state(Signal::Stokes);
        else if (tmplarch->get_npol()==4)
            tmplarch->convert_state(Signal::Stokes);
        else
            tmplarch->convert_state(Signal::Intensity);
        tmplprof = tmplarch->get_Profile(0,0,0);
    }
    catch (Error& e) {
        cerr << "Error loading template:" << endl;
        cerr << e << endl;
        return(-1);
    }

    /* Generate gaussian template if needed */
    if (gauss_width>0.0) try {
        /* Get nbin,etc from first archive */
        tmplarch = Archive::load(archives[0]);
        tmplarch->fscrunch();
        tmplarch->tscrunch();
        if (tmplarch->get_npol()==4)
            tmplarch->convert_state(Signal::Stokes);
        else
            tmplarch->convert_state(Signal::Intensity);
        tmplprof = tmplarch->get_Profile(0,0,0);
        /* Replace profile with Von Mises func */
        MEAL::ScaledVonMises vm;
        vm.set_centre(0.5);
        vm.set_height(1.0);
        vm.set_concentration(1.0/(gauss_width*gauss_width));
        for (unsigned ibin=0; ibin<tmplprof->get_nbin(); ibin++) 
            tmplprof->get_amps()[ibin] = vm.compute((double)ibin / 
                    (double)tmplprof->get_nbin());
    }
    catch (Error& e) {
        cerr << "Error creating initial template:" << endl;
        cerr << e << endl;
        return(-1);
    }

    ProfileShiftFit psf;
    smooth(tmplprof);
    psf.set_standard(tmplprof);
    if (tmplprof->get_nbin()/4 < nharm) { nharm = tmplprof->get_nbin()/4; }
    psf.set_nharm(nharm);
    Pulsar::max_harmonic = psf.get_nharm();
    //psf.set_error_method(ProfileShiftFit::MCMC_Variance);

    /* Iterate template, TOA fitting */
    unsigned it=0, nskip=0, ntot=0;
    double ttot=0.0;
    MJD t0, t1;
    bool converged=false;
    while ((it<max_it) && (!converged)) {

        /* Reset sum */
        for (unsigned ipol=0; ipol<4; ipol++) 
            sum[ipol] = new Profile(tmplprof->get_nbin());
        isum = new Profile(tmplprof->get_nbin());

        /* Reset counters */
        ntot=0;
        nskip=0;
        ttot=0.0;
        t0=99999.0;
        t1=0.0;

        /* Reset norm */
        double pnorm = 0.0;

        /* Loop over files
         * TODO: maybe we want to load them all into memory 
         * at once?  The current way might be slower but scales
         * better to very large datasets.
         */
        for (unsigned iarch=0; iarch<archives.size(); iarch++) {
            try {

                /* Load data */
                arch = Archive::load(archives[iarch]);

                /* Scrunch if needed */
                if (fscrunch) arch->fscrunch();
                if (tscrunch) arch->tscrunch();

                /* Remove profile baseline */
                arch->remove_baseline();

                /* Convert pol state */
                if (use_invint)
                    //arch->convert_state(Signal::Invariant);
                    arch->convert_state(Signal::Stokes);
                else if (arch->get_npol()==4) 
                    arch->convert_state(Signal::Stokes);
                else 
                    arch->convert_state(Signal::Intensity);

                /* Loop over subint, channel */
                for (unsigned isub=0; isub<arch->get_nsubint(); isub++) {
                    subint = arch->get_Integration(isub);
                    t0 = std::min(t0, subint->get_epoch());
                    t1 = std::max(t1, subint->get_epoch());
                    ttot += subint->get_duration();
                    for (unsigned ichan=0; ichan<subint->get_nchan(); ichan++) {

                        ntot++;

                        double shift, eshift, scale, escale, chi2;
                        prof = subint->get_Profile(0,ichan);

                        /* Skip if zero weight */
                        if (prof->get_weight()==0.0) { nskip++; continue; }

                        /* New method */
                        psf.set_Profile(prof);
                        psf.compute();
                        shift = psf.get_shift().get_value();
                        eshift = sqrt(psf.get_shift().get_variance());
                        scale = psf.get_scale().get_value();
                        escale = sqrt(psf.get_scale().get_variance());
                        chi2 = psf.get_mse();

                        /* Skip this one if snr is too low */
                        if ((psf.get_snr()<snr_cutoff) 
                                || (scale<=0.0) 
                                || (chi2<=0.0)) {
                            nskip++;
                            continue;
                        }

                        for (unsigned ipol=0; ipol<arch->get_npol(); ipol++) {

                            prof = subint->get_Profile(ipol,ichan);

                            /* Shift, scale appropriately */
                            prof->rotate_phase(shift);
                            prof->scale(scale/chi2); // Orig
                            //prof->scale(1.0/chi2);
                            //prof->scale(1.0/sqrt(chi2));

                            /* Add into sum */
                            sum[ipol]->sum(prof);

                        }
                        //pnorm += 1.0/sqrt(chi2); 
                        //pnorm += 1.0/chi2; 
                        pnorm += scale*scale/chi2; // orig
                        //pnorm += scale/chi2;
                        //pnorm += 1.0;

#if 0 
                        if (arch->get_npol()==4) {
                            pprof = subint->new_PolnProfile(ichan);
                            pprof->invint(prof);
                            isum->sum(prof);
                        }
#endif

                        /* TODO Any other stats to track?  Zap bad files? */
                    }
                }

                printf("\rIter %d Archive %d/%d Skipped %d/%d      ", 
                        it, iarch, 
                        archives.size(), nskip, ntot);
                fflush(stdout);

            }
            catch (Error &e) {
                cerr << "Error processing " << archives[iarch] << endl;
                cerr << e << endl;
                archives.erase(archives.begin()+iarch);
                iarch--;
                continue; 
            }

        }

        /* TODO: figure out convergence */

        for (unsigned ipol=0; ipol<tmplarch->get_npol(); ipol++) {

            /* Normalize */
            sum[ipol]->scale(1.0/pnorm);

            /* Replace template with sum */
            for (unsigned ibin=0; ibin<tmplprof->get_nbin(); ibin++)
              tmplarch->get_Profile(0,ipol,0)->get_amps()[ibin] 
                  = sum[ipol]->get_amps()[ibin];
        }
        isum->scale(1.0/pnorm);

        /* Smooth template */
        smooth(tmplprof);
        smooth.set_hold(true);
        for (unsigned ipol=1; ipol<tmplarch->get_npol(); ipol++) 
            smooth(tmplarch->get_Profile(0,ipol,0));
        smooth.set_hold(false);

#if 0 
        /* Normalize */
        double dc = tmplprof->sum()/(double)(tmplprof->get_nbin());
        sum->offset(-1.0 * dc);
        double pow = sqrt(sum->sumsq());
        //printf("Template power = %.6e\n", pow);
        sum->offset(dc);
        sum->scale(1.0/pow);
#endif

        /* Reset template */
        psf.set_standard(tmplprof);

        /* Increase iteration count */
        it++;
    }
    printf("\n");

    /* One more pass to determine TOAs using final template */
    if (toaoutfile!="") { 
        printf("Computing and writing final TOAs\n");
        FILE *toafile = fopen(toaoutfile.c_str(),"w");
        //FILE *snrfile = fopen("snr.dat", "w");
        for (unsigned iarch=0; iarch<archives.size(); iarch++) {
            arch = Archive::load(archives[iarch]);
            if (fscrunch) arch->fscrunch();
            if (use_invint) arch->convert_state(Signal::Invariant);
            else arch->convert_state(Signal::Intensity);
            for (unsigned isub=0; isub<arch->get_nsubint(); isub++) {
                subint = arch->get_Integration(isub);
                for (unsigned ichan=0; ichan<subint->get_nchan(); ichan++) {
                    psf.set_Profile(subint->get_Profile(0,ichan));
                    psf.compute();
                    //fprintf(snrfile, "%.6e %.6e %+.6e\n", psf.get_snr(),
                    //        psf.get_mse(), psf.get_scale().get_value());
                    if (psf.get_snr() < snr_cutoff) continue;
                    psf.toa(subint).unload(toafile);
                }
            }
        }
        fclose(toafile);
        //fclose(snrfile);
    }

    /* Replace template with sum
     * if we want orig sum out 
     */
    tmplarch->get_Integration(0)->set_epoch(0.5*(t0+t1));
    tmplarch->get_Integration(0)->set_duration(ttot);
    tmplarch->update_model();
    if (use_invint) {
        tmplarch->convert_state(Signal::Invariant);
        for (unsigned ibin=0; ibin<tmplprof->get_nbin(); ibin++)
            tmplarch->get_Profile(0,0,0)->get_amps()[ibin]
                = isum->get_amps()[ibin];
    } else {
        for (unsigned ipol=0; ipol<tmplarch->get_npol(); ipol++) 
            for (unsigned ibin=0; ibin<tmplprof->get_nbin(); ibin++)
              tmplarch->get_Profile(0,ipol,0)->get_amps()[ibin] 
                  = sum[ipol]->get_amps()[ibin];
    }

    /* Save template */
    if (tmploutfile!="") tmplarch->unload(tmploutfile);

    return(0);
}
catch (Error& e) {
    cerr << e << endl;
    return(-1);
}
