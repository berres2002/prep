FROM python:3.10

RUN apt-get update && apt-get install -y libgl1-mesa-glx

RUN pip install opencv-python-headless
COPY requirements.txt /src/requirements.txt

RUN pip install -r /src/requirements.txt

COPY setup.py /src/setup.py
COPY prep /src/prep

RUN pip install /src

CMD ["auto-prep"]