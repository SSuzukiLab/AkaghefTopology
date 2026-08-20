# 可視化トラック 引き継ぎ（2026-08-21）

宛先: Codex（実装担当）および Akaghef（起床後）。
このファイルは作業指示の正本。前提の思想は下記2文書にあり、**先に読むこと**。

---

## 0. 最初に読む順序

1. `Execution/qmd思想.md` — 研究知識基盤の縦分業。`Implementation ≠ Experiment`、
   「計算と乖離した知識は腐る」、`.qmd` は experiment artifact。
2. `topology/research/ATLAS_VISION.md` — 何を作ろうとしているか。辞書軸と観測所軸、
   値レコード、evidence-class。
3. `topology/julia/PLOT_PIPELINE.md` — 描画パイプラインの現状・撤回記録・IS一覧。

**この順序を守る理由**: 8/13–20 の作業は 1 を知らないまま進み、実験層を `.qmd` ではなく
`.jl` に置くという構造的な誤りを作った。同じ失敗を繰り返さないため、思想を先に入れる。

---

## 1. 今回のスコープ

| ID | 内容 | 状態 |
|---|---|---|
| **A** | 図を「出力」から「部品」にする（SVG 出力の整備） | 着手 |
| **B** | 検証を図の上に載せる（性質検査 + 診断表示 + 合否ボード） | 着手 |
| **C** | atlas のデータ層（値レコードのスキーマ、evidence-class 語彙） | A/B の後 |
| **E** | move の可視化（PS / MP / CP、handle slide、Turaev） | 調査から |

**D（URDiagram の可視化）は全体として対象外。** ad hoc な実装であり、今回は触らない。
これに伴い次も対象外とする:

- `src/ur_diagram.jl` および `URDiagram` 関連の Live Script
  （`C240610URD` / `C240623URD` / `URDforKnotCmpl` / `dev_URDiagram`）
- 既知のテスト失敗 `run_URDforKnotCmpl_d2_exploratory_cell`
  （`simplify2!` の後に図が `[6,-6]` へ潰れ、続く move 列が適用できない）。
  **直さない。緑にしようとしない。** スコープ外として放置する。

---

## 2. 判断が変わった点（旧文書と食い違ったら本節が優先）

### 2.1 図の受入基準

**MATLAB との一致では受け入れない。** アルゴリズムが仕様として保証する品質
（曲げ最小の直交埋め込みであり、実際に埋め込みであること）を性質検査で担保する。
最小曲げ MILP は退化しており同一目的値の別解をソルバが選ぶが、**その差異は許容**する。

具体的な規律:

- テストは**目的関数値**を比較する。bending ベクトルや座標を保存値と比較しない。
- MATLAB の図は先行実装であって正解ではない。
- `examples/` に手打ちされた `bending_numbers` は MATLAB 経路のワークアラウンドであり、
  Julia では逆に図を壊す。削除対象（§4 T2）。

### 2.2 Windows での Sage（2026-08-21 訂正）

先行文書は「Windows で MATLAB から Sage は動かない / canonical fixture は取得不能」と
書いているが、**これは誤り**。`Execution/experiment/sage/run_sage4.m`（2026-06-23）は
WSL の Sage をブリッジ経由で呼び、`plot()` の PNG を MATLAB に表示するところまで通っている。

正確な線引き:

| 主張 | 判定 |
|---|---|
| `pyenv` で WSL Sage を埋め込むのは不可能 | 成立（Linux ELF vs Windows CPython） |
| Windows から Sage を実行して図を得られない | **誤り。撤回** |
| `SageWrapper` の設計は Windows で動かない | 成立 |
| canonical fixture は取得不能 | **条件付きで誤り** |

`SageWrapper` が動かない理由は限定的で、ブリッジが `/exec` ごとに
`subprocess.run([sage, "-c", code])` で新プロセスを立てるため、`sageName` の常駐名前空間が
保てないこと。**各 exec に構築式を前置して stateless 化すれば載る。** Sage に到達するのは
`VirtualLink.m` の13メソッド。

ただし**受入基準を性質検査に変えた判断は substrate の議論から独立に成立する**ので、
2.1 は維持される。fixture が取れるようになっても、それを受入ゲートには戻さない。

---

## 3. 撤回の記録（再発防止のため残す）

同じ誤りを踏まないよう、原因と規則を明記する。

**R1 — 自己ループ図への過大一般化。** 「bending を固定すれば Sage と一致」を hopf /
borromean / whitehead の3件で確認して一般化したが、3件とも自己ループを持たなかった。
`s2xs1_calculation` では12辺中3辺が不一致。
*規則*: 構造的に別カテゴリの入力（自己ループ、非連結、仮想交差）を含めずに一致を一般化しない。

**R2 — bezier 制御点を折線として検査した。** `allocate_pos.py` は折線ではなく
**bezier 制御点**を返す（3点＝2次、4点＝3次）。制御多角形の辺で線分交差を判定しても
描かれた曲線について何も測っていない。これに基づく「Sage の欠陥を発見」は全部撤回。
**上流 Sage のバグは何も実証されていない。**
*規則*: 幾何的性質を主張する前に、座標を生成しているコードの出力契約を読む。
なお端点は曲線上に乗るので、**端点だけを見る検査（P3・次数）は fork の出力に対しても有効**。

**R3 — 供給された証拠パケットを現在の状態として扱った。** Windows のログは Feb–Mar 2026
のもので、06-23 に解決コミットが親リポジトリにあった。それを fetch せずに categorical な
不可能性を5文書へ書き込んだ。
*規則*: パケットはその日付時点の事実しか立証しない。canon に不可能性を書く前に今日までの
履歴を確認する。

---

## 4. タスク仕様

### T1 — 性質検査を `src/` に実装する【B】

新規ファイル `src/layout_invariants.jl`。**既存ファイルを編集しない**（§5 の衝突表）。

検査項目（`LinkLayout` を受けて診断リストを返す。例外を投げない）:

| ID | 内容 |
|---|---|
| P1 | 全セグメントが軸平行 |
| P3 | 各辺の折線が tail 交差点で始まり head 交差点で終わる |
| P5 | 交差点座標が相異なる |
| **P4v** | **端点の座標をグループ化したとき、どの交差点も次数ちょうど4** |
| P6 | 辺どうしが交差点以外で交わらない |
| P7 | 辺が非接続の交差点上を通過しない |

**P4v が最重要。** `tail`/`head` のメタデータを一切使わず座標だけで判定するので、
経路生成が壊れていても自分の記録で隠れられない。自己ループも自然に扱える（同じ交差点に
2端点を出すので合計4のまま）。P5 違反は次数8として、辺の消失は次数不足として現れる。
計算量も O(辺) で P6 の O(n²) より安い。

返り値は診断リスト（`PLOT_PIPELINE.md` の `Diagnostic` 相当）。severity は error / warning。

**受入**: `VLExample` の19件に対して、P3 と P4v が同じ3件
（`koda_virtual` / `s2xs1_calculation` / `koda_fig13`）を検出すること。
P1 と P5 は19件とも通ること。

### T2 — `bending_numbers` の手打ち上書きを削除する【B】

`examples/VLExample.jl:42` の `bending=[-3,0,0,0,0,0,-3,2,-2,2,-2,2]` ほか、
`examples/OgraphPlotSuite.jl` の同種の上書きを削除する。

理由: MATLAB/Sage の既定解を避けるために手で探した値であり、Julia では逆に図を壊す。
`s2xs1_calculation` は上書きを消すと P3/P4v が通る。§2.1 の基準では bending を固定する
理由がない。

**受入**: 削除後に `VLExample` が 17/19 になること（残り2件は IS7）。

### T3 — IS7 を直す【B の前提】

`src/orthogonal_layout.jl` の `_route_edges`。閉路を閉じる辺の終端が既定座標に一致するか
検証されておらず、一致しないまま沈黙する。

```
koda_virtual        edge 14 head: got (-1,-1)  want (2,2)
s2xs1_calculation   edge  8 head: got (-3,0)   want (-5,0)
koda_fig13          edge  5 head: got (1,0)    want (-1,0)
```

常に1図につき1辺、常に head 側。**`koda_fig13` が4実交差（平面化後6）で最小の再現例。**
手で全12辺を追って検算済みで、edge 5 以外は正しい。

切り分け: `_piece_lengths` が閉路辺に与える長さが悪いのか、`_route_edges` の
turn / direction の扱いが悪いのかの二択。まず後者を疑う（`crossing_rotations` の
`mod(last_direction-target_slot+3,4)` と、既発見交差点へ戻る枝で `target_slot` を
再計算していない点）。

**受入**: 3件とも P3 / P4v を通ること。かつ既に通っている16件が退行しないこと。
bending ベクトルを保存値と比較するテストを**書かないこと**（§2.1）。

### T4 — SVG 出力を部品にする【A】

`src/svg_plot.jl`。現状の問題:

- 900×600 固定で内容にフィットしない。`viewBox` を内容から決める
- 色が `#0000ff` / `#faf9f7` のベタ書き。CSS 変数に出し、ダークテーマに対応する
- `<title>` / `<desc>` が無い。クラスも無いので外からスタイルできない
- 図にメタ情報が埋まっていない（対象名・重み・代数・パラメータ）
- 交差ギャップ `0.20` が**投影前の絶対量**。図が大きいほど相対的に痩せる → 出力単位にする
- 矢印が点列の**中央インデックス**。弧長中点にする（折れ点数で位置が飛ぶ）

**受入**: 既存19図が再生成でき、T1 の検査を通ること。ライトとダーク両方で読めること。

### T5 — 診断モードの描画【B】

T1 が error を返した図について、違反箇所を図の上に描く。次数が4でない交差点を囲む、
迷子の辺を色分けする、など。合否ボード（46インスタンスの一覧）も同様。

**理由**: `VISUAL_VERIFICATION.md` の目視レビューは19件中3件を見逃している。
人間の目視より5行の端点検査の方が強かった、という事実を運用に反映する。

### T6 — 値レコードのスキーマ【C】

`ATLAS_VISION.md` §4 の `(対象, 図式, 代数, 特殊元) ⟼ 値 + evidence-class + link`。
**T1〜T5 の後に着手する。** 動くパイプラインから型を発見する。先に決めると図の都合で歪む。

evidence-class の語彙は `qmd思想.md` の `Observation → Conjecture → Claim → Proof` と
`ATLAS_VISION.md` §3 WA4 が同じことを言っているので**1つに統一する**。

### T7 — move の可視化【E】

調査から。`Manifold/VirtualLink/makeMovesData.m` と `moveData/*.mat` に既存データがある。
`.mat` は MATLAB バイナリなので Julia から読む経路を先に確認すること。
`ATLAS_VISION.md` の Area 表では move 不変性が各行の最終列に立っており、
**静的な図を並べるより研究の説明力が高い**。ただし T1〜T5 より大きいので後回しでよい。

---

## 5. 衝突表 — 触ってよいファイル

複数の担当が同時に動くため、ファイル単位で分ける。

| ファイル | 担当 | 備考 |
|---|---|---|
| `src/layout_invariants.jl` | T1（新規） | 既存に触らない |
| `src/orthogonal_layout.jl` | T3 のみ | T4 担当は触らない |
| `src/svg_plot.jl` | T4 のみ | T3 担当は触らない |
| `examples/*.jl` | T2 のみ | T3 の検証には使ってよいが編集しない |
| `src/ur_diagram.jl` | **誰も触らない** | スコープ外（D 除外） |
| `julia/*.md` | Akaghef / Claude | 実装担当は編集しない。事実の相違は報告する |
| `Manifold/`, `SW/` | **触らない** | MATLAB 側は今回対象外 |

---

## 6. 証拠の規律

このプロジェクトで実際に起きた失敗に基づく。守ること。

1. **実行していないことを「動く」と書かない。** 検証できないなら「未検証」と書く。
2. **数値を文書へ転記しない。** 生成できるものは生成する（`qmd思想.md`「計算と乖離した
   知識は腐る」）。`PLOT_PIPELINE.md` の表は現状この規律に違反しており、T6 で解消する。
3. **一般化の範囲を明示する。** 「N件で確認」と書き、確認していないカテゴリを述べる。
4. **推測で埋めない。** 埋めるときは推測と明示する（`ATLAS_VISION.md` §5 の
   「暗黙知の凍結」を再生産しない）。

---

## 7. 既知の落とし穴

- **IS4**: `set_data!(v; pd=...)` は PD ラベルを保存しない。
  `[1 2 4 3; 3 4 6 5; 5 6 2 1]` → `[5 2 6 3; 3 6 4 1; 1 4 2 5]`。
  bending や replay position を PD 行列と一緒に運ぶときは `pd_code` の返り値に索引を合わせる。
- **IS5**: 上流 Sage は MILP2 の direction バケツで modulo を最初の分岐にしか適用しない
  （`direction` が 5..7 で変数が全バケツから落ちる）。Julia は全分岐に `% 4` を適用しており、
  挙動が違う。測定範囲では `direction` は最大4で未発火。到達可能性は未証明。
- **IS6**: `planarize_virtual_gauss` が1描画につき2回走る（`orthogonal_layout` と
  `_real_arc_layouts`）。結果は一致するが二重実行。
- **planarize は未照合**: 受入対象45インスタンスのうち**22件**が
  `planarize_virtual_gauss` を通るのに、Sage との照合が1件もない。移植物で最も複雑な部分。
- **submodule ポインタ**: `AkaghefExecution` の `origin/main` が指す topology は
  `b2d32d6`(6/23) で、**Julia 移植 `00f43c0`(8/14) より前**。main から見ると `julia/` が存在
  しない。分岐はしていない（`b2d32d6` は現 HEAD の祖先）ので進めれば直る。
- **Quarto 未導入**: `.qmd` トラックの前提だが今回のスコープ（A/B/C/E）には不要。

---

## 8. 進捗（2026-08-21 夜間、Claude 実施分）

| 項目 | 状態 |
|---|---|
| Windows 記述の訂正（5文書） | 完了。R3 として `PLOT_PIPELINE.md` §2 に記録 |
| T1 性質検査 `src/layout_invariants.jl` | **完了・検証済み**。P1/P3/P4v/P5/P6/P7 |
| T2 bending 上書きの削除 | **完了・検証済み** |
| T3 IS7 | Codex 実行中 |
| T4 SVG 出力 | 実装済み・図を目視確認済み。テストスイート未実行 |
| T5 診断表示 | 未着手 |
| T6 値レコード | `scripts/check_all_layouts.jl` が JSON を吐くところまで |
| T7 move 可視化 | 未着手 |

### T2 の測定結果（削除の根拠）

**全ての手打ち上書きが有害か無意味だった。** 上書きあり／なしで性質検査を回した結果:

```
C250823MSTinv_extalg/lens_2_1           P3,P4v     → ok
C250823MSTinv_extalg/s2xs1_calculation  P3,P4v,P6  → ok
C260125MSTinv_uqsl2/lens_2_1            P3,P4v     → ok
C260125MSTinv_uqsl2/s2xs1_calculation   P3,P4v,P6  → ok
C260125MSTinv_uqsl2/s2xs1_weighted      ok         → ok   （無意味）
C250818MSTinv_uqsl2/lens_2_1            P3,P4v     → ok
C250818MSTinv_uqsl2/s2xs1_calculation   P3,P4v,P6  → ok
VLExample/s2xs1_calculation             P3,P4v,P6  → ok
```

`lens_2_1` が3箇所で壊れていたのは把握していなかった。**当初の見積り（1件）より広い。**

### 新規に判明した欠陥 IS9

性質検査を `VLExample` だけでなく全44インスタンスに回して発見。`connected_sum_2`（2群）、
`connected_sum_2_again`、`connected_sum_3` が **P6 のみ**（`connected_sum_3` は P6+P7）で
失敗する。P1/P3/P4v/P5 は通る = 接続関係は正しいのに図が埋め込みになっていない。
**IS7 とは別の欠陥。** T3 の担当は追いかけないこと。詳細は `PLOT_PIPELINE.md` IS9。

### 数値の転記をやめた

`julia/scripts/check_all_layouts.jl` が全インスタンスを検査して
`artifacts/layout_check.json` を書く。`PLOT_PIPELINE.md` §6 は件数を書かずにこれを指す。
`qmd思想.md` の「計算と乖離した知識は腐る」への対応であり、`ATLAS_VISION.md` §4 の
値レコードの種でもある。

### T4 で変えた描画規約

- 固定 900×600 の投影をやめ、**viewBox を内容の bbox に合わせた**。以降すべての長さが
  格子単位で表現でき、図の大小でギャップが痩せる問題が消える
- 交差ギャップ `0.18`、角丸半径 `0.28`、頂点 `0.115`、線幅 `0.10` — すべて格子単位
- 矢印を**弧長中点**へ（従来は点列の中央インデックスで、折れ点数で位置が飛んだ）
- `<title>` / `<desc>` / `<metadata>`、`class` 属性、CSS 変数によるダークテーマ対応。
  色は**プレゼンテーション属性と CSS クラスの両方**で出す（`rsvg-convert` は CSS 変数を
  解決しないため、属性がフォールバックになる）
- `plot_svg(::URDiagram)` と `plot_urd_structure_svg` は D 除外につき**未変更**

**既知の見た目の問題（未対処）**: 重みラベルと頂点ラベルが密な図で重なる（`koda_fig13`)。
重みが `0.0` `-1.0` と float 表示になる。どちらも T5 以降。

## 9. 未決定（Akaghef の判断待ち）

- `.qmd` の物理的な置き場所（`Execution/` 直下か `topology/experiments/` か）
- `examples/` の二重責務の割り方（実験記録と test fixture を兼ねている。
  `test/runtests.jl` が **41個の `run_*()`** を呼んでいる）
- push の可否（本引き継ぎ時点では commit のみ、push は保留）
