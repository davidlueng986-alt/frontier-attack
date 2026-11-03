# 窮舉方法：HN、JC、BSD的可終止性邊界
即使超過十億年也列方法，但不能把無界枚舉說成保證停止的雙向決定法。

- HN：按頂點數枚舉全部有限graphs（n點至多2^(n(n-1)/2)個），
  用k^n assignments判k色，再用real polynomial existential決定固定graph的
  平面unit-distance可實現性。若存在不可k色的平面有限見證，枚舉最終可找到。
  緊緻性提供該方向的有限見證；若平面可k色，這個無界負見證搜尋不保證停止。
  可枚舉有理／代數參數的週期著色作upper witness搜尋，但不聲稱涵蓋所有著色。
- JC：逐degree建立coefficients及碰撞points的complex polynomial existential系統，
  用Groebner／代數可行性決定每個固定degree。非零限制可引入逆元auxiliary。
  degree沒有已證全域界時，排除所有已跑degree不決定無界問題。
  稀疏support只是方法子族；排除子族不排除全部coefficients。
- BSD：有理Weierstrass曲線可按height枚舉；嚴格rank／analytic-order證書
  的通用可取得性不能憑空假設。有限掃描零mismatch不證全域rank equality，
  更不證full BSD。可以dovetail候選與可驗證證據，但不承諾對兩種答案都停止。
- 形式證明枚舉：可並行找P或¬P在某可枚舉axiom系統的證明；
  只對其中可證的命題有停止保證，不能假設系統對一切數學真命題完整。

完整有限問題MOLS10與CAP7的粗暴成本、安全包圍、reference程式與批次證據：
https://github.com/davidlueng986-alt/123321/blob/math-agent-mols10/docs/EXHAUSTIVE_METHODS.md

「包圍」不是尋找無限多有趣引理：每個引理應證明原完整枚舉中一大塊可安全刪去、
合併等價continuations、或解除一個根證明瓶頸；不保證一定有小型障礙。
