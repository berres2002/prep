FROM python:3.10

RUN apt-get update && apt-get install -y libgl1-mesa-glx

# RUN echo -e "import os\nlogin=os.env['LOGIN']\npassword=os.env['PASSWORD']\ntoku=os.env['TOKU']" > /src/prep/auth.py

RUN pip install opencv-python-headless
COPY requirements.txt /src/requirements.txt

RUN pip install -r /src/requirements.txt

COPY setup.py /src/setup.py
COPY prep /src/prep

RUN pip install /src

CMD ["auto-prep"]