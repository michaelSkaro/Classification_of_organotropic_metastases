FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
RUN ["apt-get", "update"]
RUN ["apt-get", "install", "--assume-yes", "--no-install-recommends", \
     "python3", "python3-pip", "default-jdk", "git"]
RUN ["git", "clone", "--branch", "package-branch", \
     "https://github.com/michaelSkaro/Classification_of_organotropic_metastases.git", \
     "/mot"]
WORKDIR /mot
RUN ["pip3", "install", "."]
RUN ["mkdir", "/demo-outputs"]
ENTRYPOINT ["python3", "-m", "mot.metastasis_pipeline", "-i", \
            "./samples/metastasis-demo/", "-o", "/demo-outputs", \
            "-w", "./lib/weka.jar", "-c", "./classes", \
            "-j", "./src/GainRatio.java"] 
