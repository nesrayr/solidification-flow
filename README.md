# solidification-flow

## Описание
Система реализует численное моделирование течения жидкости с фазовым переходом (затвердеванием) на основе
сочетания кастомного решателя OpenFOAM и теплового решения с использованием FEniCS

## Структура проекта

- app/ – клиентское Python‑приложение для моделирования фазового перехода в решателе и дальнейшего формирования КЭМ при помощи FEniCS,
- simulation/ – директория с настройками OpenFOAM‑кейса (constant, system, 0),
- solidificationFoam/ – исходники кастомного солвера OpenFOAM (C++).

## Требования

- C++ компилятор и утилиты сборки OpenFOAM (wmake или CMake)
- Python 3.10
- Conda
- Docker

## Установка и сборка

1. Сборка OpenFOAM‑солвера

```
cd solidificationFoam/solvers
wmake
```

2. Настройка окружения для FEniCS

```
conda env create -f environment.yml
conda activate fenics-env
```

Может потребоваться еще указать путь до активированной conda-среды

```
export PATH="$CONDA_PREFIX/bin:$PATH"
```

3. Запуск приложения

```
python3 client.py --json test.json --path ../simulation/constant
# в test.json заданы для теста физико-химические свойства чистой меди (Cu) в твёрдой и жидкой фазах
```

После выполнения всех шагов моделирования результаты симуляции будут лежать в simulation/VTK, полученная КЭМ модель в simulation/FEM_results
