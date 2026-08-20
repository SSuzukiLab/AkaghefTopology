# Julia topology port — 作業メモ（2026-08-14）

> ## 2026-08-20 追記：以下の記述は失効した
>
> 本文はそのまま残す（経緯の監査のため）。ただし次の点は現在の事実と異なる。
>
> - **`julia/` は追跡済み。** commit `00f43c0`（8/14 06:57）で 591 ファイルが
>   追加され、`e941b76` で origin/master とマージ済み。ワークツリーはクリーン。
> - **ランタイムは復旧済み。** 8/20 に julia 1.6.3 と node v26.0.0 を確認。
>   `Pkg.test()` は 16 pass / 1 error（`run_URDforKnotCmpl_d2_exploratory_cell`
>   の `simplify2!` 後に move 列が適用できない。これは未解決）。
> - **再開順 4 は保留（2026-08-21 訂正）。** `pyenv` で WSL の Sage を埋め込む
>   のは不可能だが、`Execution/experiment/sage/run_sage4.m`(6/23) が WSL の Sage
>   をブリッジ経由で呼んで図を返すところまで通っている。fixture は `SageWrapper`
>   を stateless 化すれば取得可能。ただし受入基準が性質検査に変わったため、
>   fixture を待つ手順自体が不要になった。詳細は `julia/PLOT_PIPELINE.md` 第7節と R3。
> - **再開順 5 の受入基準は変更された。** 図の受入は MATLAB との一致ではなく
>   性質検査（P1/P3/P5/P6）に置き換えた。最小曲げ MILP の退化に起因する
>   不一致は許容する。`julia/LIVESCRIPT_PORTS.md` の冒頭を参照。
> - **C251020 の「不合格」は撤回。** Borromean・Whitehead とも P1/P3/P5/P6 を
>   通る。MATLAB との経路差は退化した最適解の選択によるもので、現基準では
>   問題にしない。
> - **代わりに開いている課題**は IS7：`_route_edges` が閉路を閉じる辺を
>   誤った交差点で終端させる。`VLExample` の19本中3本
>   （`koda_virtual`、`s2xs1_calculation`、`koda_fig13`）で発生。
>   `s2xs1_calculation` は `examples/VLExample.jl` の bending 上書きを
>   削除すれば消える。
>
> 現在の作業状態は `julia/PLOT_PIPELINE.md` を正とする。

## 現状

- 対象は Execution 配下の `.mlx` 全60本。17本の topology/algebra 系だけを
  対象とする以前の整理は誤りで、残り43本も未移植として残っている。
- Julia 実装は `topology/julia/` にあり、Git 上は未追跡。既存の
  `Manifold/VirtualLink/VirtualLink.m` の変更は本作業では触っていない。
- 再起動後、`julia` と `node` が PATH から消えている。ランタイムを復旧・確認
  するまで、テストを通ったとは再主張しない。

## 図の受入状況

- C251020 (`projects/invariant/C251020KnotDsl2bfromUR.mlx`) の Borromean と
  Whitehead は **不合格**。MATLAB PNG と直接比較すると、外周経路、交差点の
  相対位置、矢印位置が異なる。
- Julia 自動レイアウトの構造一致は描画一致の根拠にしない。
- `julia/examples/C251020ReplayFixture.jl` と `plot_svg(...;
  replay_positions=...)` は、canonical MATLAB `REdgeTable.Position` を受けて
  arc 順、交差点端点、有限座標を検査しながら再生する。
- JSON では各 Position を N-by-2 の `[real, imag]` 行列として運ぶ。
- PNG から手計測した座標 fallback は推測であり、作成後に不適切と判断して削除。
  受入根拠に使わない。
- Windows 正本 fixture の依頼は送信済みだが、
  `Akaghef-Bridge/pending/windows_to_mac/` に C251020 応答は未着。

## 計算側の進捗

- C251030: N=4 数値還元は MATLAB 原典と同じ
  `(C,V,W)=(-1/51,[13,-13],[-113/17])`、逆トレース `-339`。
- C251030: 最終セルの恒等式
  `1/trace(D) + det(A + I - circshift(I,1)) = 0` を
  `run_C251030_final_matrix_workflow()` として追加し、再起動前に PASS。
- C251020: rank-one と multi-rank `calcInv2` を実装。有理点で exact reduction と
  記録した Laurent 式を照合している。ただし完全な多変数記号簡約は未完。
- URDforKnotCmpl: 最終 D2 の連続 move 列を
  `run_URDforKnotCmpl_d2_exploratory_cell()` に移植した。実行結果は再起動後に
  Julia ランタイム復旧をしてから確認する。`lim0` は非破壊である。

## 再開順

1. `command -v julia` と `julia --version` で Julia 1.6 系の実体を復旧・確認する。
2. Julia を**同時に1プロセスだけ**起動し、短い shard でテストする。長時間・複数
   プロセスを残さない（今回 PC が重くなった原因候補）。
3. C251030 final matrix と URDforKnotCmpl D2 entry point を個別に実行する。
4. Windows canonical `REdgeTable.Position` fixture 到着後、
   `render_c251020_replay_fixtures(json_path, output_root)` を実行し、PNG 化して
   MATLAB PNG と並べて目視確認する。差があれば renderer を修正する。
5. C251020 が合格するまで「全図確認済み」「移植完了」とは言わない。
6. 残り43 Live Script を実際に Julia 化し、各々でコード・数値・図を確認する。

## 参照

- `julia/LIVE_SCRIPT_SCOPE.md` — 60本の棚卸し
- `julia/LIVESCRIPT_PORTS.md` — セル単位の実装・未完了状況
- `julia/VISUAL_VERIFICATION.md` — 目視判定ログ
- `julia/scripts/matlab/export_replay_json.m` — canonical fixture exporter
- `julia/examples/C251020ReplayFixture.jl` — fixture importer / render
