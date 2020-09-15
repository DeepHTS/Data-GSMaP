FROM ubuntu:18.04

#ENV PYTHONUNBUFFERED 1
ENV DEBIAN_FRONTEND noniteractive
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get -y install binutils vim git python3-pip wget unzip
#RUN apt-get update -y && apt-get upgrade -y
#RUN apt-get -y install libproj-dev gdal-bin libgdal-dev
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/pip3 /usr/bin/pip
RUN pip install -U pip && hash -r pip

USER root
RUN mkdir /code
WORKDIR /code

## zsh
#RUN bash -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"
#RUN echo 'shell "/usr/bin/zsh"' >>  /home/$USER/.screenrc

# GDAL
#ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
#ENV C_INCLUDE_PATH=/usr/include/gdal
RUN apt-get update -y && apt-get upgrade -y
RUN apt install -y libproj-dev gdal-bin libgdal-dev
RUN pip install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==`gdal-config --version`

# install nodejs and npm for plotly on jupyterlab
ENV NODEJS_VERSION v12
RUN apt-get install -y curl
ENV PATH=/root/.nodebrew/current/bin:$PATH
RUN curl -L git.io/nodebrew | perl - setup && \
    nodebrew install-binary ${NODEJS_VERSION} && \
    nodebrew use ${NODEJS_VERSION}

RUN pip install -U pip
COPY requirements.txt /code/requirements.txt
RUN pip install -r /code/requirements.txt
RUN pip uninstall GDAL -y
RUN pip install GDAL==$(gdal-config --version) --global-option=build_ext --global-option="-I/usr/include/gdal"
#RUN pip install torch==1.5.0+cu101 torchvision==0.6.0+cu101 -f https://download.pytorch.org/whl/torch_stable.html
RUN pip install -U setuptools

RUN mkdir /src
RUN mkdir /dataset
WORKDIR /src
RUN chmod -R a+w .

#ENV alias python=python3.6