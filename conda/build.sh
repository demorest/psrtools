export CXXFLAGS=$(echo "$CXXFLAGS" | perl -pe 's/-std=\S+\s/-std=c++11 /')
echo "build.sh updated CXXFLAGS=${CXXFLAGS}"
make && make install
