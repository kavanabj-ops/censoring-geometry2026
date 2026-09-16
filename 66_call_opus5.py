# -*- coding: utf-8 -*-
# 66_call_opus5.py — 调 opus5 确认 A1b 结果（加 thinking disabled 防截断）
import json, urllib.request, time

API_KEY = "sk-nQTBFYWHECeIMIfZ7522C8A70b9c4d179bCb2eD2253fFc5d"
BASE = "https://api.apiyi.com/v1"
OUTDIR = r"D:\paper\02_analysis\进程2对话"

CONTEXT = """你是进程2（审稿专家 claude-opus-5）。进程1 跑完了 A1b 描述性普查（GDSC2 真实数据）。

【A1b 结果】
- 靶向药（228个，PUTATIVE_TARGET 为基因符号）：sufficient 711 组合，X%（WT格删失率>0.90）=37.4%（>0.80 为 47.1%）
- 细胞毒药（20个，机制词如 DNA crosslinker/microtubule）：sufficient 4209，X%=21.2%（>0.80 为 32.2%）
- insufficient（逐格 n<5）：靶向 94.4%，细胞毒 82.0%
- 删失不平衡 r_wt-r_mut：中位≈0（两药都如此）

【进程1 诚实自评】
1. X% 靶向(37%)>细胞毒(21%)，差 1.76 倍，但没达到你预期的"远高于"。
2. insufficient 94.4%（靶向）可能是比 X% 更硬的发现："逐格 n 不足"比"删失"更普遍。
3. 删失不平衡≈0 是我用"非沉默突变"（含非热点）定义 mut 格导致的伪影——BRAF 只有 V600E 敏感，G596R 等非沉默突变不敏感，混入后拉高 mut 格删失率。你 round3 要求"只纳 OncoKB Oncogenic 或 hotspot/truncating"，我偷懒用了非沉默突变。

【请简短回答，每条≤150字，直接给结论，不要 extended thinking、不要长推导】
1. 突变定义要不要收紧到 hotspot？CCLE 无现成 OncoKB 注释时，最快怎么识别 hotspot？还是 X% 只看 WT 格、非沉默突变已够用、不必重跑？
2. X% 37% vs 21% 和 insufficient 94% 这两个数字在论文里怎么定位？insufficient 94% 是否升级为新卖点？
3. 主卖点还是"可操作规则 R1-R5"吗？"""

payload = {
    "model": "claude-opus-5",
    "messages": [{"role": "user", "content": CONTEXT}],
    "max_tokens": 16000,
    "temperature": 0.2,
    "thinking": {"type": "disabled"},
    "stream": True,
}
req = urllib.request.Request(BASE + "/chat/completions",
    data=json.dumps(payload).encode(),
    headers={"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"})

t0 = time.time()
print("开始请求 ...")
resp = urllib.request.urlopen(req, timeout=600)
full = ""
n_chunks = 0
last_t = time.time()
for raw in resp:
    line = raw.decode("utf-8", errors="replace").strip()
    if not line.startswith("data:"):
        continue
    data = line[5:].strip()
    if data == "[DONE]":
        break
    try:
        chunk = json.loads(data)
        delta = chunk["choices"][0].get("delta", {})
        if "content" in delta and delta["content"]:
            full += delta["content"]
            n_chunks += 1
            last_t = time.time()
        # 也检查是否有 reasoning/thinking 字段（诊断用）
        if "reasoning" in delta or "thinking" in delta:
            pass
    except Exception:
        pass
    if time.time() - last_t > 60:
        print(f"  [进度] 60秒无新 content，已收 {len(full)} 字符")
        last_t = time.time()

elapsed = time.time() - t0
print(f"耗时 {elapsed:.1f}s，返回 {len(full)} 字符（{n_chunks} 个 content chunk）")
if full:
    with open(OUTDIR + r"\round6_a1b.md", "w", encoding="utf-8") as f:
        f.write(full)
    print("=" * 70)
    print(full)
else:
    print("空返回——thinking 可能仍在占用 token，需进一步禁用")
