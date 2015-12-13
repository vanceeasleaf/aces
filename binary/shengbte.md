tar xvf spglib.tar.gz
./configure --prefix=/a/lib/spglib/
make install
cp /a/lib/splig/lib/libsymspg.so /a/lib/splig/lib/libsymspg.so.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/a/lib/splig/lib/
tar xvf shengbte.tar.gz
cp arch.make.example Src/arch.make.example
edit spglib path
make 

