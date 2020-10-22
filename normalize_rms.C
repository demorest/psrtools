/***************************************************************************
 *
 *   Copyright (C) 2008 by Paul Demorest
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

using namespace std;

#include "Pulsar/Application.h"
#include "Pulsar/StandardOptions.h"
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/ProfileShiftFit.h"
#include "Pulsar/PhaseWeight.h"
#include "strutil.h"
#include "math.h"

using namespace Pulsar;

//
//! Archive rms normalizer
//
class normalize_rms : public Pulsar::Application
{
public:

    //! Default constructor
    normalize_rms ();

    //! Setup
    void setup ();

    //! Process the given archive
    void process (Pulsar::Archive*);

protected:
    void add_options (CommandLine::Menu&);

    bool do_weight;
};


normalize_rms::normalize_rms ()
  : Application ("normalize_rms", "Scale data so offpulse RMS = 1.0")
{
    do_weight = false;
}

void normalize_rms::add_options (CommandLine::Menu& menu)
{
    CommandLine::Argument* arg;
    arg = menu.add (do_weight, 'w');
    arg->set_help("apply scaling to weights rather than to data");
}



void normalize_rms::setup ()
{
}

void normalize_rms::process (Pulsar::Archive* archive)
{

    // Output filename
    cout << "# " << archive->get_filename() << endl;

    Reference::To<Archive> copy = archive->clone();
    copy->tscrunch();
    copy->convert_state(Signal::Intensity);

    bool scale_indep_subints = true;

    // Get baseline values
    // Use this for applying single scaling vs freq to all subints
    vector< vector< Estimate<double> > > mean;
    vector< vector< double> > var;
    if (scale_indep_subints == false)
        copy->get_Integration(0)->baseline_stats(&mean, &var);

    // Loop over subints
    for (unsigned isub=0; isub<archive->get_nsubint(); isub++) {

        // Use this to get a new scaling for each subint
        if (scale_indep_subints)
            archive->get_Integration(isub)->baseline_stats(&mean, &var);

        // Loop over channels
        for (unsigned ichan=0; ichan<archive->get_nchan(); ichan++) {

          double varnorm = fabs(sqrt(var[0][ichan])/mean[0][ichan].val);
          if (mean[0][ichan].val==0.0) varnorm=0.0;

          if (do_weight) {
            double weight;
            if (varnorm>0.0) {
                weight = 1.0 / sqrt(var[0][ichan]);
                if ((archive->get_state()==Signal::Coherence) 
                        || (archive->get_state()==Signal::PPQQ))
                    weight = 1.0 / sqrt(0.5*(var[0][ichan]+var[1][ichan]));
            } else 
                weight = 0.0;
            if (archive->get_Integration(isub)->get_weight(ichan)!=0.0) 
                archive->get_Integration(isub)->set_weight(ichan, weight);
          } else {
            // Scale orig data
            double scale;
            if (var[0][ichan]>0.0) {
                scale = 1.0 / sqrt(var[0][ichan]); // Single-pol or stokes
                if ((archive->get_state()==Signal::Coherence) 
                        || (archive->get_state()==Signal::PPQQ))
                    scale = 1.0 / sqrt(0.5*(var[0][ichan]+var[1][ichan]));
                for (unsigned ipol=0; ipol<archive->get_npol(); ipol++) 
                {
                    //scale = 1.0 / sqrt(var[ipol][ichan]); // Separate pols
                    archive->get_Profile(isub,ipol,ichan)->scale(scale);
                }
            } else {
                archive->get_Integration(isub)->set_weight(ichan, 0.0);
            }

            // Reset all non-zero weights to 1
            // XXX why did I do this.... ?
#if 0 
            if (archive->get_Integration(isub)->get_weight(ichan)!=0.0) 
                archive->get_Integration(isub)->set_weight(ichan,1.0);
#endif

          }
        }
    }

    // Unload corrected archive
    // TODO use standard unload options
    archive->unload(replace_extension(archive->get_filename(), "norm"));

}

static normalize_rms program;

int main (int argc, char** argv)
{
    return program.main (argc, argv);
}

