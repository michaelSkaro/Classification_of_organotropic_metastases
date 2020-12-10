FROM ubuntu
COPY  ./ /mot
ENV DEBIAN_FRONTEND=noninteractive
RUN ["apt-get", "update"]
RUN ["apt-get", "install", "--assume-yes", "--no-install-recommends", \
     "python3", "python3-pip", "default-jdk", "git"]
#RUN ["git", "clone", \
#     "https://github.com/michaelSkaro/Classification_of_organotropic_metastases/tree/package-branch", \
#     "/mot"]
WORKDIR /mot
RUN ["pip3", "install", "-r", "requirements.txt"]
ENTRYPOINT ["/bin/bash"]
