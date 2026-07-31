"""관측된 |군간 차이| 중앙값이 '잡음만 있을 때의 기대치'를 넘는지 본다.

효과가 0 이어도 |차이| 의 중앙값은 0 이 아니다. 표본 잡음 때문에 양수가 나온다.
그 기대치는 표준오차로 정해진다.
  신뢰구간 폭 = 2 * t_crit * SE   ->  SE = 폭 / (2 * t_crit)
  효과가 0 일 때 |차이| 의 중앙값 = 0.6745 * SE   (반정규분포의 중앙값)
이 값을 '잡음 바닥' 으로 삼아 관측값과 비교한다.
"""
import pandas as pd, numpy as np, os
from scipy import stats

ROOT = r'C:\Users\louis\Desktop\kaist\26 여름\KASA'
R = pd.read_csv(os.path.join(ROOT, 'results/tables/stage3d_osd21_clock.csv'))
R = R[R.gene.notna()].copy()

LAB = {"HLU+Reloaded vs NormalLoaded": "HLU+재하중 vs 정상하중",
       "Flight vs GroundControl":      "우주비행 vs 지상대조",
       "HLU vs NormalLoaded":          "HLU vs 정상하중",
       "HLU+Reloaded vs HLU":          "HLU+재하중 vs HLU"}

print(f"{'대비':<24}{'n':>7}{'관측':>8}{'잡음바닥':>9}{'배수':>7}{'CI폭':>7}{'0배제':>7}")
print("-" * 70)
out = []
for k, lab in LAB.items():
    s = R[R.comparison == k].copy()
    na, nb = int(s.n_a.iloc[0]), int(s.n_b.iloc[0])
    df = na + nb - 2
    tc = stats.t.ppf(0.975, df)
    se = (s.ci_high - s.ci_low) / (2 * tc)
    floor = float(np.median(0.6745 * se))       # 효과 0 일 때 |차이| 중앙값
    obs = float(s.diff_log2.abs().median())
    nsig = int(((s.ci_low > 0) | (s.ci_high < 0)).sum())
    w = float((s.ci_high - s.ci_low).median())
    out.append((lab, f"{na}v{nb}", obs, floor, obs / floor, w, nsig))
    print(f"{lab:<24}{na}v{nb:<5}{obs:>8.3f}{floor:>9.3f}{obs/floor:>7.2f}x{w:>7.2f}{nsig:>5}/15")

print("\n해석")
print("  '배수' 가 1.0 이면 관측된 차이가 잡음만으로 기대되는 크기와 같다는 뜻이다.")
print("  즉 그 대비에서는 실제 효과의 증거가 없다.")

o = dict((r[0], r) for r in out)
rl = o["HLU+재하중 vs HLU"]; fl = o["우주비행 vs 지상대조"]
print(f"\n  재하중 3.5시간 대비: 관측 {rl[2]:.3f} / 잡음바닥 {rl[3]:.3f} = {rl[4]:.2f}배")
print(f"    -> 잡음 기대치와 사실상 같다. 3.5시간이 변화를 '일으켰다'고 말할 수 없다.")
print(f"  우주비행 대비:       관측 {fl[2]:.3f} / 잡음바닥 {fl[3]:.3f} = {fl[4]:.2f}배")

print("\n  '배수' 순위와 0배제 개수 순위 비교")
for lab, n, obs, fl_, r, w, ns in sorted(out, key=lambda x: -x[4]):
    print(f"    {lab:<24} 배수 {r:.2f}x   0배제 {ns}/15")
