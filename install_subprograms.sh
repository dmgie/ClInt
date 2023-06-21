# Here we either build them from source (i.e minimap2) if cloned as a submodule
# or we download them from the web (i.e. samtools, bcftools, etc.) as a binary

# minimap2
# check if in programs/minimap2 or if minimap2 is in the path
if [ ! -x ./programs/minimap2/minimap2 ] || command -v minimap2 >/dev/null 2>&1; then
    git clone --depth=1 www.github.com/lh3/minimap2 programs/minimap2
    cd programs/minimap2
    make
    cd ../..
fi

# fastp
if [ ! -d ./programs/fastp ] || command -v fastp >/dev/null 2>&1; then
    # Latest linux binary
    wget http://opengene.org/fastp/fastp -O programs/fastp
    chmod a+x ./fastp
fi

# samtools
if [ ! -d ./programs/samtools-1.17 ] || command -v samtools >/dev/null 2>&1; then
    wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 -O programs/samtools.tar.bz2
    tar -xjf programs/samtools.tar.bz2 -C programs/
    rm programs/samtools.tar.bz2
    cd programs/samtools-1.17
    ./configure --prefix=../ # To install into programs folder
    make
    make install
    cd ../..
fi

# bcftools
# This is needed for samtools
if [ ! -d ./programs/bcftools-1.17 ] || command -v bcftools >/dev/null 2>&1; then
    wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 -O programs/bcftools.tar.bz2
    tar -xjf programs/bcftools.tar.bz2 -C programs/
    rm programs/bcftools.tar.bz2
    cd programs/bcftools-1.17
    ./configure --prefix=../ # To install into programs folder
    make
    make install
    cd ../..
fi

# htslib
# This is needed for bcftools/samtools
if [ ! -d ./programs/htslib-1.17 ] || command -v htsfile >/dev/null 2>&1; then
    wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 -O programs/htslib.tar.bz2
    tar -xjf programs/htslib.tar.bz2 -C programs/
    rm programs/htslib.tar.bz2
    cd programs/htslib-1.17
    ./configure --prefix=../ # To install into programs folder
    make
    make install
    cd ../..
fi



