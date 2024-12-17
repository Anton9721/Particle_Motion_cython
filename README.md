Сборка проекта с нуля:
1. Удалить все файлы, кроме .vscode, include, sourse, CmakeLists.txt, main.cpp, particle_motion.pyx, Cython_test.ipunb, setup.py

2. Проверить, есть ли в окружении рабочая версия python и cython.

3. В строке bash (или powershell) прописать команды:

mkdir build

cd build
cmake ..

cmake --build .

cd..

python setup.py build_ext-- inplace
