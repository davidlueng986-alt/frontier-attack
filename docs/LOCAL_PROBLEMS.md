# Frontier 局部關口：只作主問題閉合的手段
2026-11-03。現有新舊報告不在本次全面重新認證；不聲稱任何新數學結果。

| ID | 要回答的精確問題 | 足以收斂甚麼 | 驗收與界限 |
|---|---|---|---|
| F-H0 | 所引用的平面色數baseline真正成立且適用嗎？ | 判定6/7是否完整目標集合 | 原文statement與完整證明依賴；不靠README |
| F-H1 | 指定abstract gadget是否真的強迫different而非equal？ | 關閉錯誤抽象構造族 | 可著色正控制、端點同色假設的UNSAT證書 |
| F-H2 | 指定graph有精確平面unit-distance embedding嗎？ | 分開抽象色數與幾何障礙 | exact algebraic coordinates、互異性、每邊d²=1；非浮點 |
| F-H3 | gadget composition是否推出主目標色數？ | 防止只重證較弱下界 | 完整推導、替換接口與幾何相容性；明列>=6或>=7 |
| F-H4 | 候選全域著色含所有點與邊界嗎？ | upper-bound proof所缺介面 | 區域聯集／重疊／所有距離1情況完整覆蓋 |
| F-J0 | counterexample=>normalized system的必要方向正確嗎？ | 確定emptiness可排除的整類 | field、degree端點、support、全部corners的原文審核 |
| F-J1 | inverse transform真的回到polynomial ring嗎？ | 決定反例構造路是否可用 | 所有負冪消失；清分母後重驗det與碰撞 |
| F-J2 | 指定完整系統有精確emptiness證書嗎？ | 可關閉該完整類或提高degree bound | 例如Σh_i f_i=1獨立展開；vertices-only不能推全文support |
| F-J3 | 所有degree是否有完整覆蓋／一般證明？ | 無界JC的必要根義務 | GAP；有限D不能自動提升 |
| F-B0 | 同一curve的algebraic rank上下界可否相等？ | rank證書瓶頸 | 嚴格independence及descent/Selmer upper；不只database |
| F-B1 | 低階L導數真為零而第k階嚴格非零嗎？ | analytic-order證書瓶頸 | interval含零不是等於零；誤差與定理假設完整 |
| F-B2 | 是否只研究rank equality還是full BSD？ | 防止根命題被偷換 | full BSD額外leading coefficient義務另列GAP |

反例優先：先嘗試推翻足以終止大量支線的universal shortcut。
若成立，再研究可重用障礙與安全剪枝；不重複堆積相近local證書。
本次不開啟三題研究或任何automation；最新分支歷史保留。
