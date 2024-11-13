# research
コンパイル方法
```
/research/source$ cmake -S . -B build
/research/source$ cmake --build build
/research/source/build$ ./main "モデル名.obj"
```

新しいファイルを追加して，それを新たにコンパイルに加えたい場合は，
`/source/CMakeLists.txt`の`add_executable`内にファイル名を追加してください．