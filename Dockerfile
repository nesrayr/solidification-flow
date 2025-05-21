# Используем dev-образ с инструментами сборки
FROM opencfd/openfoam-dev:2412

# Копируем residuals.cfg
RUN mkdir -p /usr/lib/openfoam/openfoam2412/etc/caseDicts/postProcessing/numerical/
COPY simulation/residuals.cfg /usr/lib/openfoam/openfoam2412/etc/caseDicts/postProcessing/numerical/residuals.cfg
RUN chmod 644 /usr/lib/openfoam/openfoam2412/etc/caseDicts/postProcessing/numerical/residuals.cfg

# Установить переменные окружения
ENV WM_PROJECT_DIR=/usr/lib/openfoam/openfoam2412
ENV PATH=$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/bin:$PATH

# Установить сборочные тулзы
RUN apt update --fix-missing \
  && apt install -y \
  python3-venv \
  python3-pip \
  build-essential cmake flex bison git nano \
  calculix-ccx calculix-cgx gfortran \
  && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/ccx2pv \
  && /opt/ccx2pv/bin/pip install --upgrade pip \
  && /opt/ccx2pv/bin/pip install vtk ccx2paraview

# Добавляем скрипты из venv в PATH
ENV PATH=/opt/ccx2pv/bin:$PATH

# Копируем решатель и симуляцию
COPY solidificationFoam /home/openfoam/solidification-flow/solidificationFoam
COPY simulation /home/openfoam/solidification-flow/simulation

RUN mkdir -p /home/openfoam/OpenFOAM/user-v2412/platforms/linuxARM64GccDPInt32Opt/bin && \
  mkdir -p /home/openfoam/OpenFOAM/user-v2412/platforms/linuxARM64GccDPInt32Opt/lib

RUN printf 'source $WM_PROJECT_DIR/etc/bashrc\n' > /etc/profile.d/openfoam.sh

ENV FOAM_USER_DIR=/home/openfoam/OpenFOAM/user-v2412
ENV FOAM_USER_APPBIN=$FOAM_USER_DIR/platforms/linuxARM64GccDPInt32Opt/bin
ENV FOAM_USER_LIBBIN=$FOAM_USER_DIR/platforms/linuxARM64GccDPInt32Opt/lib

# Задаём точку входа
CMD ["/bin/bash"]
