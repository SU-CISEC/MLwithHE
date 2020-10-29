# CLion remote docker environment (How to build docker container, run and stop it)
#
# Build and run:
#   docker build -t ml_with_he:1.0 -f Dockerfile .
#   docker run -d --cap-add sys_ptrace -p127.0.0.1:2222:22 --name ml_with_he ml_with_he:1.0
#   ssh-keygen -f "$HOME/.ssh/known_hosts" -R "[localhost]:2222"
#
# stop:
#   docker stop clion_remote_env
#
# ssh credentials (test user):
#   user@password

FROM ubuntu:latest

RUN DEBIAN_FRONTEND="noninteractive" apt-get update \
  && apt-get install -y ssh \
      build-essential \
      autoconf \
      cmake \
      git \
      gcc \
      g++ \
      gdb \
      rsync \
      tar \
      python \
  && apt-get clean

RUN cd /

RUN ( \
    echo 'LogLevel DEBUG2'; \
    echo 'PermitRootLogin yes'; \
    echo 'PasswordAuthentication yes'; \
    echo 'Subsystem sftp /usr/lib/openssh/sftp-server'; \
  ) > /etc/ssh/sshd_config_test_clion \
  && mkdir /run/sshd

RUN useradd -m user \
  && yes password | passwd user

CMD ["/usr/sbin/sshd", "-D", "-e", "-f", "/etc/ssh/sshd_config_test_clion"]