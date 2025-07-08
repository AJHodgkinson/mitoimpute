# Use CentOS as the base image
FROM centos:latest

RUN cd /etc/yum.repos.d/
RUN sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-*
RUN sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*

#Load dependencies

RUN dnf install -y redhat-rpm-config
RUN yum -y update
RUN yum -y --enablerepo=extras install epel-release

RUN yum install -y \
git \
python2 \
wget \
epel-release \
python2-pip \
gcc \ 
python2-devel \
make \
zlib-devel \
gcc-c++ \
bzip2 \
bzip2-devel \
ncurses-devel \
xz-devel \
perl-Env \
java-devel \
perl-core \
gsl \
gsl-devel

RUN yum -y install dnf-plugins-core && \
    yum config-manager --set-enabled powertools

RUN yum -y install R

#Get plink
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
RUN unzip plink_linux_x86_64_20230116.zip
RUN mv plink /bin/

RUN wget -O- https://cpanmin.us | perl - App::cpanminus
RUN cpanm --mirror-only --mirror https://cpan.metacpan.org/ Statistics::R

RUN yum -y install curl-devel
RUN yum -y install libxml2-devel
RUN yum -y install openssl-devel

ENV R_LIBS_USER /usr/lib64/R/library
RUN echo 'options(repos = list(CRAN = "https://cloud.r-project.org/"))' > /usr/lib64/R/etc/Rprofile.site

RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20231212.zip
RUN unzip plink2_linux_x86_64_20231212.zip
RUN mv plink2 /bin/

RUN wget https://github.com/gusevlab/fusion_twas/archive/master.zip
RUN unzip master.zip
RUN cd fusion_twas-master
RUN mv master.zip old.zip

RUN wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
RUN tar xjvf LDREF.tar.bz2

RUN wget https://github.com/gabraham/plink2R/archive/master.zip
RUN unzip master.zip

RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('RColorBrewer')"
RUN Rscript -e "install.packages('glmnet')"
RUN Rscript -e "install.packages('methods')"

RUN wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
RUN gunzip gemma-0.98.5-linux-static-AMD64.gz
RUN mv gemma-0.98.5-linux-static-AMD64 gemma
RUN mv gemma fusion_twas-master/
RUN chmod +x fusion_twas-master/gemma

RUN Rscript -e "install.packages('plink2R-master/plink2R/',repos=NULL)"

#Get liftover
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
RUN mv liftOver /bin/
RUN chmod +x /bin/liftOver
RUN wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
RUN wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz


RUN ln -s ./ output
