# Lagrangean 粒子追跡シミュレータの開発

名前はまだない. 仮で"PARTRAS" (Particle TRAcking Simulator)と名付けている. 

候補:
- DROPS (DROplet Particle Simulator)
- LAPTRAN (Lagrangean Approach for Particle Tracking in FortRAN)
- FoPAS (Fortran-based Particle Advection SImulator)


## やりたいこと
- 流れ場を浮遊する粒子の運動を計算するためのツールを提供する.
- 拡張性を考慮して, クラスを用いた設計を行う.
- 他の言語に対するインターフェースを提供する. 

## fortran coding rule
- ソースファイルは UpperCamelCase で命名する
- fortranの変数や関数, モジュールなどは原則 snake_case とする.
- 変数名については原則省略形を用いない 
  - 慣例的な語句は省略形を用いても良いが, コメントで何の意味かを明記する.
- コメントは FORD スタイルに従う.
- 派生型の名前は`_t`を末尾につける
- モジュールの名前は`_m`を末尾につける
- private変数は`_`を末尾につける
- モジュール変数は`mv_`を先頭につける. 
- 後始末は `delete_`を先頭につける. 
- クラスのコンストラクタは関数ではなくサブルーチンとする. 
- コンストラクタは `construct_`を先頭につける. 
 
### do OOP or not ?
OOPにすると, 
- 😄 継承が使えるので再利用が楽. 
- 😄 インタープリターでメンバの候補が出てくるので楽.
- 😢 メンバの隠蔽がめんどくさい. いちいちゲッタがいる. 
  - 隠蔽しなくても可能. Pythonの思想に近い.
  - Fortranの場合subroutineのintent属性である程度制限できる. 冗長になりがち. 
- 😢 そもそも, シミュレータの場合, OOPの恩恵が少なく, 自己満足になりがち.
  - 継承とかあまり使わないことが多い. MDみたいに複数の積分法や様々な統計処理が関わって大規模かつ複雑化するなら恩恵はあるが, 
    飛沫計算の場合はやるべき処理がほぼ決まっているので, むしろ無駄に複雑化する恐れがある. 

そこで...

シミュレーションで重要な2つのデータ: 粒子情報と流れ場情報 `particle_data`, `flow_field`はモジュール変数とする. 
* とくに `flow_field`はprotected属性をつけ, モジュール外からはメンバ変数は参照onlyとする. メンバの変更はモジュールに定義した関数のみ許可. 
  グローバル変数みたいなものだが, protected属性で外部からは参照しかできない構造になっていることに注意. 
* `particle_data`は書き込み, 参照が自由. **紳士協定**として, `particle_data`の初期化, 破棄, データ追加はモジュール定義サブルーチンから行う. 
  * クラスにすればええやんと思うかもしれないが, 構造体にすることでPythonとのやり取りが簡単になる利点がある. 
* いわばモジュールそのものが一つのクラスのインスタンスのようなイメージ. 

その他
* 動的リスト`Idlist`, 格子読み込み処理`ugrid_importer_t`など, 一時オブジェクトとして利用される可能性の高いもの, 別の方法を試したい可能性のある処理はクラスにする. ほかに updater, dump_writerなど. 

## シミュレータの中身
**主なコンポーネント**
- `ugrid_importer_t`
  - 何らかのフォーマットからデータを読み込んで`ugrid_struct_t`を構築する. 
  - read_file: ファイルを読み込む. `ugrid_struct_t`を返す (subroutine)
  - delete: データを消す. 
- `ugrid_struct_t`
  - 外部ファイルから読み取った格子データをメモリに保持するための構造体. 
  - `flow_field_t`の構築はこの構造体を引数に渡して行う. こうすることで, コンストラクタを一つに統一できる. ファイルフォーマット依存性を消すための処置. 
- `particle_data_t`
  - 粒子の位置, 速度
  - 生死状態
  - 書き出しと読み取りの定義
  - 粒子の追加
- `flow_field_t`
  - 節点座標, セル中心, 面中心, 境界面, 隣接関係
  - 流体の速度場
  - セル速度 --> 節点速度への変換 (必要?)
- `motion_t`
  - 粒子の運動の定義
- `simulator_t`
  - シミュレーション全体. 
    - `particle_data_t`と`flow_field_t`はすでに構築済みとする.
  - メインループ:
    - 流れ場データの更新 
    - 粒子の運動
    - コールバック
    - データの書き出し
  - コールバックの登録, 呼び出し
- `field_updater_t`
  - 流れ場を更新するためのヘルパー.
  - 指定した`importer_t`で更新する. 
  - 粒子運動の時間刻み幅と流れ場の時間刻み幅から, アップデートのタイミングを計算する  
- `callback_t`
  - `simulator_t`で呼び出すコールバック処理. 

**データ書き出し**
- `flow_field_t`, `particle_data_t`の生データを出力する. 
- タイムステップごとに別のファイルで出力するか, HDF5みたいな圧縮フォーマットで, シミュレーションのデータをすべて
  ひとつのファイルに書き出すか (HOOMDでいう`.gsd`みたいなやつ.)

可視化ファイルは Legacy VTKとする. 当分困ることはないのでこれで. ただのモジュールサブルーチンで十分. 

**データファイルのフォーマット**
計算結果の出力ファイルのフォーマット. 拡張子は`.pdata`とする. 
- タイムステップごとにファイルを出力する. ファイル名は`tag_#####.pdata`とする. 
- `tag`はユーザーが指定する. 数字`#####`は整数で, シミュレーションの時間ステップとする. 

```markdown
# number of particle
(number of particle)
# particle attributes (pos, vel, f, ...)
(data for 1st particle)
(data for 2nd particle)
(...)
```

- バイナリ形式の場合, 文字列は書き込まない. 
  - number of particle $N_{p}$ (4byte)
  - particles (56byte*$N_p$) ※現在の定義

- 書き込み形式は`stream`で.
- 可視化ファイルの作成は別にする. 

**リスタート処理**
- `particle_data_t`: ファイル`.pdata`から取得.


誰がやる? --> DumpReader?
- DumpReaderはオブジェクトを初期化できる特別な存在になっちゃう.
- 後処理で使うこと考えると別に問題ないか.

処理の流れ: 
1. `backup.pdata`を読み込む. 
2. `backup.pdata`の粒子情報から`pdata`を初期化する
3. `backup.pdata`のファイル名から流れ場を構築する
4. 時間ステップ`ncyc`は誰が持つ?
updater, writerなどのオブジェクトはユーザが逐一用意する. 

**CFDメッシュサポート**
- `ugrid_importer_t`を継承して作成する. 
- `read_file`を各自実装する. 
  - `.pre`, `.fld`, `.fph`, CUBE関連のやつ
- テスト用にLegacy VTKのインポータを用意している. 

**外部モジュールへの依存**
- stdlib: ロガーやエラーハンドル, システム関連は使えるかも. 
- VTKFortran: 不要. 
- HDF5: 有力候補

*廃止予定*
- `struct_array3_t` (deprecated)
  - ベクトルx,y,z成分をメンバに持つ構造体
  - 座標データは構造体の配列で実装する.
  - ４つ目の成分をフィールドに加えるなどしても良い. 
  - 作っては見たが, やっぱり使わない. 使いにくいので.

## motion
`particle_data_t`にフィールドを追加せず, こちらで例えば粒子の半径や粒子に働く力を設定する. 
```Fortran

type motion_base_t
  
  contains
  procedure integrate_a_step
  procedure update_status

end type

```

`integrate_a_step`で　粒子の位置と速度を更新する. 独自のフィールドがある場合はこれを更新することもできる. 

```Fortran
pure subroutine integrate_a_step(this, r, v, u, dt)
  class(motion_base_t),intent(in) :: this
  type(scalar4_t),intent(inout) :: r
  type(scalar4_t),intent(inout) :: v
  type(scalar4_t),intent(in   ) :: u
  real(dp),intent(in) :: dt

end subroutine

```

`update_status`はかなり自由に粒子の状態を編集できる. 強制的に粒子を止めることもできる. 粒子半径の変化もここでやる. 
この処理の妥当性の保証は, このクラスを拡張したユーザーに責任がある. 

## 追加の処理
例えば, 粒子間の合体処理, 粒子の追加処理. 
これは `simulator_t` のコールバックで処理する. 

## src/
ウイルス飛沫シミュレータ. 処理をモジュールごとにわけて書いておく. メイン関数で呼び出す

- system.f90: ウイルス飛沫シミュレータ本体. main関数. 
  - mv_field_updaterの呼び出し
  - mv_dump_writerの呼び出し
- init.f90: 流れ場(mv_flow_field)の初期化, 粒子データ(mv_pdata)の生成
  - 粒子データの生成: 与えた中心点周りにランダムで点を打つ. 速度はゼロ. 半径は共通. 
- update.f90: 流れ場の更新
- config.f90: コンフィグファイルの読み込みと自動生成
- dump.f90: ダンプライターの初期化

## Pythonとの連携
pythonのインターフェースを提供する. 
- python側ではラッパークラスを呼ぶ
- ラッパークラス内部ではfortranのインターフェース関数を呼び出す. 
- fortranのインターフェース関数内部ではクラスのインスタンスを確保する.
- particle_dataの中身はすべてfortran側で確保する. 

ユーザ定義のcomputeクラスなどはどうやってインスタンスを作成する?
考えているのは以下. `system`クラスに`add_compute`を用意し, `compute`クラスのポインタを保持させる.
ひょっとしたら, このルーチンはいらないかも. Pythonに公開する用のサブルーチン内で`add_compute`を呼び出すようにすれば済む. 

なお, `add_compute`では, `compute`クラスとその継承クラスのポインタを保持する. ポインタの配列はFortranの規格では配列に対するポインタとなってしまうので, 構造体で代用する. 

```Fortran
type compute_holder
    class(compute),pointer :: ptr_ => null()
end type
```
これをメンバとして保持する

```Fortran
type system
    type(compute_holder),allocatable :: computes(:)
    !...
end type

```

Pythonでは

```Python
class compute:
  def __init__(self, args*):
    # ...
    # ctypeslibか何かでfortran共有ライブラリを呼び出す
    fort_lib.init_compute(...)
    pass

  def hoge(self, xx):
    # 必要であれば他のメソッド
    fort_lib.set_params(...)
    pass
  
```

FortranではPython側で作ったインスタンスに対応するミラークラスのインスタンスが作成される

```Fortran
module compute_something_m
use system, only : global_system !systemクラスインスタンスのグローバル変数

!computeクラスの定義
type,extends(compute) :: compute_something
! 諸々の定義

end type

type(compute_something) compute_somthing_mv

contains
subroutine init_compute(...) bind("C", name="init_compute")
  ! 必要な引数を羅列

  ! 諸々の処理
  call compute_somthing_mv%init(...)

  ! systemにこのクラスのインスタンスを追加
  call global_system%add_compute(compute_something_mv)
end subroutine

! ...
end module

```

面倒なのは, `compute`の派生クラスを作成するたびにこれを作る必要があること.
あと, `allocatable`属性が含まれるクラスは, インスタンス作成の場所に注意. 
サブルーチン内部で作ったクラスは`allocatable`属性の変数は自動的に解放されるため.
モジュール変数にしないといけない.

motionも同様にやればいい．