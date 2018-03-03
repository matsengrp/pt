FROM debian:stretch

RUN apt-get update && \
    apt-get install -y build-essential

RUN apt-get install -y autoconf \
                       bison \
                       cmake \
                       file \
                       flex \
                       libtool

COPY . /pt
WORKDIR /pt
RUN make && make test

RUN cp _build/src/ptw_threads /usr/local/bin
