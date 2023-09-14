FROM ubuntu:20.04

RUN apt-get update -y

# manually enable gcc8 installation 
RUN echo 'deb [arch=amd64] http://archive.ubuntu.com/ubuntu focal main universe' >> /etc/apt/sources.list

# make sure gcc, python3 (dev) and pip are up to date
RUN apt-get update && apt-get install -y \
    gcc-8 g++-8 python3 python3-pip \ 
    python3-dev 

RUN pip install mnnpy nxfvars papermill scanpy

# downgrade package 
RUN pip install scipy==1.8.1

CMD ["python3"]


