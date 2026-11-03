# 唯一成功條件：主問題 OPEN → CLOSED
2026-11-03。本次只更新規約與研究入口，沒有重跑或關閉本repo三題。

局部定理、可重現程式、負結果、文獻量、checker passes和分支數都是工具。
唯一成功是相應主問題被完整證明／否證，並經獨立驗證。
不能把局部問題已解、候選方法已否定，說成主問題成功或完成百分比。

## 精確根命題
- HN：決定平面色數K，需要χ>=K與χ<=K。採6/7二擇一前先核實該baseline；
  一個兩距離圖不是unit-distance見證；數值距離近似不是精確幾何。
- 平面JC：任意C[x,y] Keller map是否可逆；反例可用真正多項式、常數非零Jacobian
  和兩個不同點的相同像。有限degree結果是子目標，不是無界JC。
- BSD：本repo現有rank mismatch目標是rank equality部分，不能冒稱full BSD。
  真正反例需同一E/Q的嚴格代數rank及analytic order，且兩者不同。

現有checker的degree/support、數值精度或backend限制是實作範圍，不是數學真假的定義。
若有效見證超出工具能力，應擴充驗證器；不能把工具拒絕當成數學反證。

## 局部問題如何服務根
每項工作列parent、精確量詞、必要／充分方向、完整coverage、input版本、
證據、獨立checker、負控制、budget和PASS/FAIL/UNKNOWN後續。
圖以AND表示缺一不可、OR表示替代充分方法；GAP不得假裝atomic task。
局部排除只能刪其精確範圍，方法失敗不等於反命題成立。
保留完整暴力／半判定基線，無論成本是否超十億年，但清楚標示不保證終止的方向。

入口：[局部問題](LOCAL_PROBLEMS.md)、[枚舉方法邊界](EXHAUSTIVE_METHODS.md)。
跨專案已識別24項清單及實作：
https://github.com/davidlueng986-alt/123321/blob/math-agent-mols10/docs/LOCAL_PROBLEM_CATALOG.md
