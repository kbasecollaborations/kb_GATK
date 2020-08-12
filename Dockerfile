FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get -y install wget 
RUN mkdir -p /kb/module/deps
WORKDIR /kb/module/deps
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
RUN unzip gatk-4.1.3.0.zip

RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.0/picard.jar
RUN apt-get update && apt-get install -y bwa \
    samtools \
    tabix \
    r-base


RUN R -e "install.packages('gplots', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"



# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
