# ParTraS - Particle Trajectory Simulator

![Fortran Standard](https://img.shields.io/badge/Fortran-2008%2F2018-blue)
<img src="https://img.shields.io/badge/-Python-F9DC3E.svg?logo=python&style=flat">
![License](https://img.shields.io/github/license/hirotori/partras)

ParTraS は、非構造格子で計算された流れ場を浮遊する微粒子や飛沫の追跡を目的とした Lagrange 粒子シミュレーションフレームワークです。
高精度な計算と柔軟な拡張性を提供し、研究や工学分野での使用を想定しています。

---

## 主な機能 (Features)

- 非構造格子に対応
- Legacy VTK formatの格子データの入力が可能
- ほかの非構造格子データもインポート可能
- 飛沫, エアロゾルのシミュレーションに対応
- 100% Modern Fortranで記述
- Pythonインターフェース `PyPTS` によるシミュレーションの実行

### PyPTS

`PyPTS`はParTraSのPythonインターフェースです. Pythonから簡単にシミュレーションが実行できるように, `ParTraS`の主要なライブラリへのインターフェースを提供します. シミュレーションは呼び出されたFortranライブラリによって行われるので, 純粋なPythonライブラリに比べて高い計算パフォーマンスが期待されます. 

`PyPTS`には以下の機能が用意されています. 
- 粒子データの初期化, 取得
- 流れ場の初期化および更新
- シミュレーション実行
- 飛沫, ウイルス飛沫, エアロゾルに適した運動方程式


---

## インストール方法 (Installation)

### 必要な環境
- Fortran コンパイラ: 以下のコンパイラで動作が確認されています
  - gfortran
  - ifort
- Git
- CMake

### 手順
1. リポジトリをクローンします:
   ```bash
   git clone https://github.com/hirotori/partras.git
   cd partras
   ```
2. ビルドディレクトリを生成し, cmakeコマンドでプロジェクトを構成する. 
   ```bash
   mkdir build
   cd build
   cmake ..
   ```
   オプションでビルド方法とコンパイラを設定できる. 
   - `-DCMAKE_BUILD_TYPE`: `debug`, `release`
   - `-DCMAKE_Fortran_COMPILER`: `gfortan`, `ifort`

    `gfortran`を使って`release`でビルドする場合
    ```bash
   cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_Fortran_COMPILER=gfortran
    
    ```


3. `make`コマンドでコンパイルを実行する
   ```bash
   make
   ```

4. `PyPTS`を使用する場合は, 環境変数にパスを追加してください. 
   
   ```bash
   export PYTHONPATH=/Users/username/partras/":$PYTHONPATH"
   ```