# -*- coding: utf-8 -*-
# 80_call_astra.py - 进程3 独立审稿（gpt-6-astra），评估投 3-4 分杂志录用可能性
import json, urllib.request, time, os

API_KEY = "sk-nQTBFYWHECeIMIfZ7522C8A70b9c4d179bCb2eD2253fFc5d"
BASE = "https://api.apiyi.com/v1"
OUTDIR = r"D:\paper\02_analysis\进程3对话"
os.makedirs(OUTDIR, exist_ok=True)

with open(r"D:\paper\02_analysis\manuscript\00_full_draft.md", encoding="utf-8") as f:
    draft = f.read()

CONTEXT = """你是独立审稿人（进程3，GPT-6-Astra），与进程2（claude-opus-5）完全独立、不受其结论影响。

【背景】进程1 完成了一篇方法学论文《Censoring-geometry collapse of mutation x lineage interaction estimands in pharmacogenomic panels》。进程2 已审 12 轮，认为可投 Briefings in Bioinformatics（IF ~6.8），预期 major revision 后接受。

【你的任务】换一个角度独立评估：这篇稿子【降档投 3-4 分杂志】的录用可能性。请逐条回答：

1. 录用概率：投 3-4 分杂志（如 BMC Bioinformatics IF~3.0、BMC Genomics IF~4.0、Database IF~3.0、Briefings in Functional Genomics 等），录用概率给一个百分比区间，并说明理由。

2. 主要风险点：哪些地方可能被 3-4 分杂志的审稿人拒掉？逐条列出（越具体越好，指出稿件里的具体段落/缺陷）。

3. 若投 3-4 分：是"可以直接投"还是"仍需补什么"？需要补的话，补什么最划算？

4. 对比：与 6-8 分（CSBJ 6.0 / BiB 6.8）相比，3-4 分是否明显更稳妥？录用概率差多少？

5. 明确建议：投哪个档、哪个具体杂志、预期结局（接收 / major revision / reject）。

请直接给结论，每条简洁（≤200字），不空话、不客套。

【稿子全文】
""" + draft

payload = {
    "model": "gpt-6-astra",
    "messages": [{"role": "user", "content": CONTEXT}],
    "max_tokens": 16000,
    "temperature": 0.2,
    "stream": True,
}
req = urllib.request.Request(BASE + "/chat/completions",
    data=json.dumps(payload).encode(),
    headers={"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"})

t0 = time.time()
print("开始请求 gpt-6-astra ...")
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
    except Exception:
        pass
    if time.time() - last_t > 60:
        print(f"  [进度] 60秒无新 content，已收 {len(full)} 字符")
        last_t = time.time()

elapsed = time.time() - t0
print(f"耗时 {elapsed:.1f}s，返回 {len(full)} 字符（{n_chunks} chunk）")
if full:
    with open(OUTDIR + r"\round1_astra_3-4分.md", "w", encoding="utf-8") as f:
        f.write(full)
    print("=" * 70)
    print(full)
else:
    print("空返回——可能需要禁用 thinking 或调整参数")
