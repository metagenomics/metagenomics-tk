FROM python:3.12.0a7-windowsservercore-ltsc2022
ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
RUN apt update && apt install -y procps
