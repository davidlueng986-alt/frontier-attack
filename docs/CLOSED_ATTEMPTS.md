# Closed attempts：唯一終點是主問題 OPEN → CLOSED

2026-11-03。三題均未被本專案解決。局部方法被排除、程式通過或找到較弱結果，
都不是主目標成功。本輪不設定成功機率或完成百分比，也不啟用任何 automation。

## 第一個動作：使終端判定真的對應主命題

現有 `tools/accept.py` 是 legacy gate，不能單獨支持 solved：

- HN 只檢查 proof path 存在，沒有實際重播 DRAT/LRAT；點座標也必須明確限於平面。
- JC 要先確認兩個真正的 polynomial components、兩個 variables、精確有限點；
  不應靠未核實的 degree 125 或某個支撐形狀來定義反例。
- BSD 的 method/check 字串不是 rank 或 analytic-order 證書；
  不可把未證的 algebraic rank parity 當成排除有效證據的規則。

先補這些實際驗證接口與負控制，再接受根結論。這是防止假閉合的必要工作，
不把修好 gate 當成三題進展的替代品。

## HN：直接取得根見證，或完整上界證明

根問題是決定 χ(R²)。如果已核實 6≤χ≤7，才可收窄為：

1. 下界路：顯式有限平面 unit-distance graph，不可六色，附完整可重播 UNSAT trace。
   配合既有全平面七色上界即 exact 7。只得到不可五色，仍未閉合。
2. 上界路：定義全平面六色著色，完整覆蓋所有點、重疊及邊界，
   並證明所有距離一點對異色。配合已驗六色下界才得 exact 6。

repo 引用的 OpenAI `Five.lean` entrypoint 已讀到，的確有 unrestricted no-five-coloring
statement；本環境尚未完成 Lean build 或 dependency axiom audit，不能只由檔案存在
把下界依賴標成獨立重驗完成。也不能直接把其中 Fin 5 換成 Fin 6 當新證明。

### 下一個直接構造 attempt

- 只生成真正有精確代數座標的平面 graph；不以兩距離 graph 代替。
- SAT 的主問題始終是六色。SAT 回傳一份可行著色，將其作為下一批結構生成
  要排除的反例；幾何候選來自經核對的代數變換與圓交點，而不是容差邊。
- 單一新頂點不可能從一個五色著色直接強迫七色；需要成批、互相制約的結構。
- 若能證明整個候選生成族總可六色，就停止該族。既有只達四色的glue不再放大。
- 任一不可六色候選立即停止擴搜尋，轉入 exact geometry、獨立編碼與 trace check。

這是 heuristic candidate discovery，不宣稱搜尋完整，也不預先承諾可找到見證。
另一條全域證明路要逐項列出從 no-six-coloring 到 exact 7 的完整推導；
六色週期構造失敗只排除該 ansatz，不否定全平面六色。

## 平面 Jacobian：真正 polynomial 碰撞，或全 degree 的證明

反例的充分終點是 F=(P,Q)∈C[x,y]²、det JF=c≠0、p≠q、F(p)=F(q)，全部精確。
不要求預設 degree pair，不把 Laurent／rational pair 當真正 polynomial map。

任何碰撞反例都可經可逆 affine domain/target transformations 化成：

    G(0,0)=G(1,0)=(0,0), det JG=1。

理由：令 domain affine map 把 0,e1 送到 p,q；其 determinant δ≠0。
再平移 target 並用固定線性變換把 cδ 正規化為1。這保持 polynomial 性與 max degree。
所以每個固定 degree D 有一個完整有限係數方程問題；逐 D 搜尋可作反例半判定，
沒有已證 degree 上界時，有限 D 的 emptiness 不是全域JC。

### 已執行的第一個 direct attempt

從 repo 可精確驗算的三維 collision source 出發，對列出的三對碰撞點，
完整檢查所有通過該兩點的複數 affine planes 和所有 constant linear output projections。
多項式鏈式法則與 minors 的 exact Bezout 證書排除了整個構造族。
見 [精確報告](../findings/2026-11-03-jc-affine-projection.md)。
這不是平面JC證明；它確切關閉了這個降維方法，避免再做有限參數取樣。

### 下一個 direct attempt

允許 collision-preserving polynomial embeddings，例如對每個座標：

    H_i(u,v)=p_i+(q_i-p_i)u+a_i(u²-u)+b_i v+c_i uv+d_i v²。

它保證 H(0,0)=p、H(1,0)=q。組合已知三維 source 再作兩個線性輸出，
得到的二維 map 自動是 polynomial 且保留碰撞；剩下直接解 det J=1 的完整係數系統。
degree≤2 embedding 的12個係數加3個 output-minor weights 是有限代數問題。
若找到候選，須重建實際兩個 linear forms，並獨立驗算二維 determinant、碰撞與座標。
若有完整 emptiness 證書，只排除這個 bounded embedding family；timeout 是UNKNOWN。
不再以 vertices-only／leading-edge 系統可解，推論 full system 可解或可逆重建。

全域 proof 路則必須給所有 degree 的推導橋樑或有效歸約界；沒有此橋樑，
再多 (72,108)/(75,125) 局部排除都不稱為全域證明。

## BSD：用嚴格分離界線直接否證，不強求兩邊 exact

full BSD 包含 rank equality 及 leading coefficient 等算術內容。
證明 rank equality 不是證明 full BSD；一個真正 rank mismatch 則足以否證 BSD。

更短的反例終點是以下任一組**已驗證界線**：

    rank E(Q) ≥ rL > kU ≥ ord_{s=1}L(E,s)，或
    ord_{s=1}L(E,s) ≥ kL > rU ≥ rank E(Q)。

因此不必同時求兩個 exact integers。令 Λ(s)=h(s)L(E,s)，h在1附近解析且非零，
則兩者 vanishing order 相同。嚴格證明 Λ^(j)(1)≠0 已給 analytic order≤j；
若又有 rank≥j+1 的獨立證據，就是足以閉合的反例，不必再證低階導數全零。

### 下一個 direct attempt

- 選有可驗證獨立有理點／rank lower certificate 的曲線，不單靠database欄位。
- 只計算足以與該 rank bound 分離的導數；functional-equation 強制為零的階直接跳過。
- 對有限和、尾項、積分與 rounding 全給嚴格誤差；數值近零和 precision rerun 都不算證書。
- 若得到不重疊 bounds，立即轉向獨立 arithmetic backend 與 analytic checker。
- 若 intervals 始終重疊，報 UNKNOWN／未找到反例；不把擴大 conductor sweep 當一般證明。

環境目前有gcc，但gp/sage/mwrank尚未安裝，舊Windows工具限制不是永久數學障礙。
先確認能取得何種嚴格 lower／upper 證據，再選曲線family；不先花算力做無法認證的掃描。
若走「BSD為真」路，必須有適用所有 E/Q 的新一般性證明；有限曲線驗證不能替代。

## 每輪交付

只回報：精確命題、完整量詞／範圍、原始見證或證書、獨立驗證、
根問題是否真的閉合，以及剩餘哪個阻塞。三題保持分開，不用一題的局部結果替另一題交差。
