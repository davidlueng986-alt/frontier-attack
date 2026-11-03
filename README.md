# frontier-attack

最新執行計畫：[Closed attempts](docs/CLOSED_ATTEMPTS.md)。三題主問題仍OPEN。
首輪已精確排除一個直接Jacobian降維構造族，見 [證據與範圍](findings/2026-11-03-jc-affine-projection.md)。

**唯一成功條件：相應主問題 OPEN → CLOSED。** 局部定理與工具都是手段，不能代替 solved。

研究入口：[成功規約](docs/MISSION.md) · [局部關口與驗收](docs/LOCAL_PROBLEMS.md) · [完整／半判定窮舉邊界](docs/EXHAUSTIVE_METHODS.md)。

以下既有記錄保留；其外部baseline、成功機率與checker限制不等於已驗根證明。

公開 repo：**https://github.com/davidlueng986-alt/frontier-attack**

三個 2026 年仍然開放的具體問題，加上一套**反替換**的 agent 工作流。

這不是「用 10,000 個 agent 硬解千禧年難題」的專案。這是一個**檢查器優先**的
工作台：先把「什麼算解決」寫成程式，再讓 agent 去產生能被程式接受的物件。

> 沒有 checker 的 agent，就是在寫散文。

## 三個目標

| # | 問題 | 精確目標 | 檢查器 | 目前狀態 |
|---|---|---|---|---|
| 1 | Hadwiger–Nelson | 決定 χ(ℝ²) = 6 還是 7 | `hn.py`（精確算術 + 精確著色數） | 6 ≤ χ ≤ 7；已知最小 5-chromatic 圖 509 點 |
| 2 | 平面 Jacobian 猜想 | n = 2 時 det JF = const ⇒ F 可逆 | `jac.py`（符號驗證 + 碰撞搜尋） | n ≥ 3 已被否定（2026-07）；n = 2 開放 |
| 3 | BSD 猜想 | 找 rank ≠ analytic rank 的 E/ℚ | `bsd.py`（不變量守門員） | rank ≤ 1 已知；一般情形開放 |

這三個問題的難度**不是同一種**：

- **#1 是搜尋問題。** 只要寫出一個 7-chromatic 的有限 unit-distance graph，就等於
  證明 χ(ℝ²) = 7（de Bruijn–Erdős 緊緻性）。物件可被機器在毫秒內驗證。
- **#2 是半搜尋問題。** 反例只要寫下來就能驗證；但 n = 2 的已知結構讓反例的
  存在空間非常窄。這是三者中最難用算力推進的。
- **#3 是資料問題。** 反例是一條曲線加兩個整數，驗證需要獨立後端；但真正的
  瓶頸是 rank 的計算，不是搜尋。

## 快速開始

```bash
python tools/verify.py                     # 跑完所有自我檢查
python problems/hadwiger_nelson/hn.py      # 精確 Moser spindle: n=7 m=11 chi=4
python problems/hadwiger_nelson/search.py --mode glue --trials 100
python problems/jacobian/jac.py            # 驗證器自我測試
python problems/jacobian/enumerate.py      # 清一個有限箱子
python problems/bsd/bsd.py                 # 守門員自我測試
```

## 目前跑出來的結果（可重現）

- Moser spindle：7 點、11 邊、χ = 4，**精確符號驗證通過**。
- Glue 搜尋（30 次、每次最多 60 點）：找到 34 點 / 69 邊 / χ = 4 的圖，
  120 位精度重新驗證通過。沒有找到 χ = 5。
- 二次 Jacobian 映射枚舉：15,625 個係數組合中 13 個 det 為常數，
  **0 個反例** —— 與 Wang (1980) 的定理一致。
- 稀疏三次族：625 個組合中只有 1 個 det 為常數（恆等映射），
  說明「稀疏多項式」這個搜尋空間會塌縮成三角映射。

這些都是**負結果**，而且是我刻意要的：它們證明檢查器不會亂放行。

## 目錄

```
docs/WORKFLOW.md        think → research → build → verify → scale（先想，再查，再做）
docs/ORCHESTRATION.md   替換現象、角色、admissibility 閘門、tangent tax
docs/SCORECARD.md       賠率、kill criteria、什麼算成功
plan/*.md               Phase 0 產物：terminal object / 拆解 / 合格動作
research/*.md           Phase 1 產物：物件與數字，附來源
problems/*/             三個問題的檢查器與搜尋器
agents/*.md             THINKER / RESEARCHER / GENERATOR / BREAKER / CARTOGRAPHER
findings/*.md           每個 pilot 的結果與「為什麼這個空間不行」
runs/                   每次執行的 JSON 產物
tools/verify.py         一鍵跑完所有檢查
```

## 工作流（強制順序）

1. **THINK** — 一頁：terminal object、terminal check、牆、3–5 個 load-bearing
   子問題、合格動作清單。清單為空就寫 `NO_ADMISSIBLE_ACTION` 並停。
2. **RESEARCH** — 背景並行，每人一個問題，回傳物件與數字，不寫程式。
3. **BUILD** — 每個 agent 一個產物，成功 = 終端 checker 回 0。
4. **VERIFY** — 只用兩個：A 跑 checker，B 用**不同方法**重寫一次。不一致就停。
   不要第三個。
5. **SCALE** — 只有在量到 `p > 0` 之後，才放大到 `N ≈ 1/p`。

三個 pilot 的 `p` 都是 0。所以那三條分支不該放大，該換表示法。

## 這個專案的核心主張

你觀察到的「AI 一直去解無關的小引理」不是模型缺陷，是**規格缺陷**：

1. 「研究 X」沒有成功謂詞，agent 只好用「產出可驗證的東西」當代理指標。
2. 引理局部可驗，目標不可驗 —— 理性的 agent 會優化可驗的那個。
3. 長上下文會讓目標衰減。
4. 沒有對抗者去殺掉離題分支。
5. 最貴的一個：**用數學結構（引理）而不是用產物（候選物件）來切分工作**。

所以本專案的規則只有一條：**每個 agent 擁有一個物件，成功是 checker 回傳 0。**

細節見 `docs/ORCHESTRATION.md`。

## 不做什麼

- 不聲稱解決任何一個問題。
- 不把「Lean 編譯通過」等同於「數學社群接受」。
- 不把數值搜尋的 negative result 說成定理。
- 不在沒有獨立後端的情況下接受 BSD 的 rank 數字。
