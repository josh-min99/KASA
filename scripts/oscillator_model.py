"""
Week 3-4: 중력 입력을 가진 분자시계 진동자 모델.

모델
  Goodwin/Gonze 계열 코어 음성 되먹임 루프 3변수.
  NR1D1(REV-ERBα) 단백질이 BMAL1 전사를 억제한다.

    B : BMAL1 활성      m : Nr1d1 mRNA      P : NR1D1 단백질

    dB/dt = v1 * K1^n/(K1^n + P^n) - v2*B/(K2 + B) + u(t)
    dm/dt = k3*B                   - v4*m/(K4 + m)
    dP/dt = k5*m                   - kdeg*P/(K6 + P)

  중력이 두 경로로 들어간다 (선행연구 대응):
    경로 1 (전정계, 급성·위상의존) : u(t) — 중력 자극 펄스의 일시적 전사 입력
        Fuller 2020 (doi:10.1038/s41598-020-65496-x): 2G 펄스가 체온 리듬 재동조를
        가속하고, 양측 전정계 손상 쥐에서는 효과가 사라진다.
    경로 2 (미토파지, 만성·진폭의존): kdeg — NR1D1 분해율
        IJMS 2024 (doi:10.3390/ijms25094853): 모사 미세중력에서 SCN 미토파지가 결핍되고
        NR1D1 '단백질'만 증가(mRNA 불변) → 단백질 분해율 저하로 모델링한다.

  두 경로를 동시에 자유롭게 두면 식별 불가능하므로 시간 스케일로 분리해 고정한다.

파라미터 주의
  Hill 계수 n=8 은 Goodwin 계열의 알려진 제약이다(3변수 루프가 안정 한계순환을 가지려면
  높은 협동성이 필요). 다량체 형성·다단계 인산화를 하나로 뭉뚱그린 값으로 해석해야 하며,
  정량적 예측이 아니라 정성적 구조 예측(PRC 모양, 동조 영역의 존재)에만 쓴다.
  파라미터는 20-28h 주기 격자 탐색에서 상대진폭이 최대인 조건으로 선택했다.

산출: data/model_*.csv
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data")
os.makedirs(OUT, exist_ok=True)

PAR = dict(v1=0.7, K1=1.0, n=8.0, v2=0.50, K2=1.0,
           k3=0.7, v4=0.35, K4=1.0,
           k5=0.7, K6=1.0)
KDEG0 = 0.35
DT = 0.02
Y0 = np.array([0.5, 0.6, 0.7])

# 위상 기준 변수. 중력 펄스 u(t) 는 dB/dt 에 직접 더해지므로 B 에는 펄스 자체가
# 인공 peak 를 만든다. 그 B 로 위상을 재면 세기가 커질수록 주기 추정이 붕괴한다
# (초기 구현에서 Arnold tongue 의 주기비 1.0 부근이 비는 형태로 드러났다).
# 하류 변수인 NR1D1 단백질(P, 인덱스 2)로 위상을 잰다. 웻랩 관측 대상과도 일치한다.
PHASE_VAR = 2
# 마우스 SCN 에서 Nr1d1 은 대략 CT6 에 정점.
CT_AT_PHASE_PEAK = 6.0


def rhs(y, u, kdeg, p=PAR):
    B, m, P = y[0], y[1], y[2]
    dB = p["v1"] * p["K1"] ** p["n"] / (p["K1"] ** p["n"] + P ** p["n"]) - p["v2"] * B / (p["K2"] + B) + u
    dm = p["k3"] * B - p["v4"] * m / (p["K4"] + m)
    dP = p["k5"] * m - kdeg * P / (p["K6"] + P)
    return np.stack([dB, dm, dP])


def integrate(y, kdeg, n_steps, u_fn=None, t0=0.0, record=None):
    """벡터화 RK4. y shape (3, N). record: None | 'B' | 'full'."""
    if record == "B":
        out = np.empty((n_steps, y.shape[1]), dtype=np.float32)
    elif record == "full":
        out = np.empty((n_steps, 3, y.shape[1]), dtype=np.float32)
    else:
        out = None
    t = t0
    for i in range(n_steps):
        u1 = u_fn(t) if u_fn else 0.0
        u2 = u_fn(t + DT / 2) if u_fn else 0.0
        u3 = u_fn(t + DT) if u_fn else 0.0
        k1 = rhs(y, u1, kdeg)
        k2 = rhs(y + DT / 2 * k1, u2, kdeg)
        k3_ = rhs(y + DT / 2 * k2, u2, kdeg)
        k4 = rhs(y + DT * k3_, u3, kdeg)
        y = np.maximum(y + DT / 6 * (k1 + 2 * k2 + 2 * k3_ + k4), 1e-9)
        if record == "B":
            out[i] = y[0]
        elif record == "full":
            out[i] = y
        t += DT
    return (y, out) if record else y


def peaks(trace_1d, t0=0.0):
    """포물선 보간 국소 최대 시각."""
    a, b, c = trace_1d[:-2], trace_1d[1:-1], trace_1d[2:]
    idx = np.where((b > a) & (b > c))[0]
    if len(idx) == 0:
        return np.array([])
    d = a[idx] - 2 * b[idx] + c[idx]
    sh = np.where(np.abs(d) > 1e-12, 0.5 * (a[idx] - c[idx]) / np.where(np.abs(d) > 1e-12, d, 1.0), 0.0)
    return t0 + (idx + 1 + sh) * DT


def to_limit_cycle(kdeg_arr, burn=1000.0):
    N = len(kdeg_arr)
    y = np.repeat(Y0[:, None], N, axis=1)
    return integrate(y, kdeg_arr, int(burn / DT))


# ------------------------------------------------------------- 1) 기저 진동
def characterize(kdeg_arr, span=500.0):
    """주기 / B진폭 / 상대진폭 / NR1D1 주기평균.

    NR1D1 은 적분 종료 시점 값이 아니라 관측 구간 평균을 써야 한다.
    종점 값은 그 순간의 위상에 따라 달라져 kdeg 효과를 가린다.
    """
    y = to_limit_cycle(kdeg_arr)
    y, tr = integrate(y, kdeg_arr, int(span / DT), record="full")
    N = len(kdeg_arr)
    per = np.full(N, np.nan); amp = np.full(N, np.nan)
    rel = np.full(N, np.nan); pm = np.full(N, np.nan); pk_ = np.full(N, np.nan)
    for j in range(N):
        b = tr[:, 0, j]
        p = tr[:, 2, j]
        pks = peaks(b)
        if len(pks) >= 4:
            per[j] = np.median(np.diff(pks[1:]))
        amp[j] = (b.max() - b.min()) / 2
        rel[j] = amp[j] / max(b.mean(), 1e-9)
        pm[j] = p.mean()          # NR1D1 주기 평균
        pk_[j] = p.max()          # NR1D1 최대
    return per, amp, rel, pm, pk_


# ---------------------------------------------------- 2) kdeg 스윕 (미세중력)
def kdeg_sweep(tau0):
    kd = np.round(np.linspace(0.20, 0.50, 31), 4)
    per, amp, rel, pm, pmax = characterize(kd)
    df = pd.DataFrame({"kdeg": kd, "period_h": per, "amplitude": amp,
                       "rel_amplitude": rel, "NR1D1_mean": pm, "NR1D1_max": pmax})
    df["oscillating"] = df.rel_amplitude > 0.05
    print("\n[2] NR1D1 분해율(kdeg) 스윕 — 미세중력 = 미토파지 결핍 = kdeg 저하")
    print(df.iloc[::3].round(4).to_string(index=False))
    base = df.iloc[(df.kdeg - KDEG0).abs().argmin()]
    for frac in [0.9, 0.8, 0.7, 0.6]:
        r = df.iloc[(df.kdeg - KDEG0 * frac).abs().argmin()]
        print(f"    kdeg {int(frac*100)}%: 주기 {r.period_h:6.2f} h  "
              f"상대진폭 {r.rel_amplitude:6.3f} ({(r.rel_amplitude/base.rel_amplitude-1)*100:+6.1f}%)  "
              f"NR1D1평균 {r.NR1D1_mean:.3f} ({(r.NR1D1_mean/base.NR1D1_mean-1)*100:+.1f}%)  "
              f"{'진동' if r.oscillating else '진동소실'}")
    dead = df[~df.oscillating]
    if len(dead):
        print(f"    → 진동 소실 임계: kdeg < {dead.kdeg.max():.3f} "
              f"({dead.kdeg.max()/KDEG0*100:.0f}% of baseline)")
    print("    IJMS 2024 관찰(NR1D1 단백질 축적 + 리듬 진폭 감소)과 방향이 일치")
    df.to_csv(os.path.join(OUT, "model_kdeg_sweep.csv"), index=False, encoding="utf-8-sig")
    return df


# ------------------------------------------------------------------- 3) PRC
def compute_prc(tau, kdeg=KDEG0, strength=0.30, dur=1.0, n_phase=48, horizon_cycles=40):
    """한계순환 위상별 펄스 인가 → 점근 위상 이동."""
    y1 = to_limit_cycle(np.array([kdeg]))
    # 한 주기 분량의 상태 궤적을 기록해 위상별 초기조건을 만든다
    _, full = integrate(y1, np.array([kdeg]), int(2 * tau / DT), record="full")
    bt = full[:, PHASE_VAR, 0]
    pk = peaks(bt)
    if len(pk) < 2:
        return None
    i0 = int(round(pk[0] / DT))
    step = (pk[1] - pk[0]) / n_phase
    idxs = np.clip((i0 + np.arange(n_phase) * step / DT).astype(int), 0, full.shape[0] - 1)
    Y = full[idxs, :, 0].T.astype(float).copy()      # (3, n_phase)

    kd = np.full(n_phase, kdeg)
    horizon = int((dur + horizon_cycles * tau) / DT)

    def u_fn(t):
        return np.where(t < dur, strength, 0.0)

    _, fp = integrate(Y.copy(), kd, horizon, u_fn=u_fn, record="full")
    _, fc = integrate(Y.copy(), kd, horizon, record="full")
    trp, trc = fp[:, PHASE_VAR, :], fc[:, PHASE_VAR, :]

    # 위상 이동은 ±tau/2 부근에서 감김(wrapping) 모호성이 생긴다.
    # 펄스 직후부터 주기별 이동량을 추적해 연속적으로 이어붙인 뒤 점근값을 취한다.
    shifts = np.full(n_phase, np.nan)
    converged = np.zeros(n_phase, bool)
    for j in range(n_phase):
        pp, pc = peaks(trp[:, j]), peaks(trc[:, j])
        k = min(len(pp), len(pc))
        if k < 8:
            continue
        tau_j = np.median(np.diff(pc[1:]))
        d = pp[:k] - pc[:k]                     # 주기별 원시 이동량
        # 첫 주기 이동량을 (-tau/2, tau/2] 로 정규화한 뒤 연속 unwrap
        d[0] = ((d[0] + tau_j / 2) % tau_j) - tau_j / 2
        # 인접 차이가 tau/2 를 넘으면 tau 만큼 보정
        for i in range(1, len(d)):
            while d[i] - d[i - 1] > tau_j / 2:
                d[i] -= tau_j
            while d[i] - d[i - 1] < -tau_j / 2:
                d[i] += tau_j
        tail = d[-4:]
        shifts[j] = float(tail.mean())
        converged[j] = bool(np.ptp(tail) < 0.05 * tau_j)   # 점근 수렴 여부

    phase_h = np.arange(n_phase) * tau / n_phase
    return pd.DataFrame({
        "internal_phase_h": np.round(phase_h, 3),
        "CT_h": np.round((phase_h / tau * 24 + CT_AT_PHASE_PEAK) % 24, 2),
        "phase_shift_h": shifts,
        "converged": converged,
    })


# ------------------------------------------------------------ 4) Arnold tongue
def arnold_tongue(tau, kdeg=KDEG0, dur=1.0, n_s=21, n_t=31):
    # 생리적으로 의미 있는 약~중간 세기 구간을 촘촘히 본다.
    # 더 강한 자극에서는 고차(n:m) 잠금과 진폭 사멸이 섞여 1:1 동조 판정이 무의미해진다.
    strengths = np.linspace(0.0, 0.20, n_s)
    periods = np.linspace(tau * 0.82, tau * 1.18, n_t)
    Sg, Tg = np.meshgrid(strengths, periods, indexing="ij")
    S, T = Sg.ravel(), Tg.ravel()
    N = len(S)
    kd = np.full(N, kdeg)
    y = to_limit_cycle(kd, burn=600.0)

    def u_fn(t):
        return np.where((t % T) < dur, S, 0.0)

    T0 = 150 * tau
    y = integrate(y, kd, int(T0 / DT), u_fn=u_fn)                    # 과도상태 소거
    _, tr = integrate(y, kd, int(60 * tau / DT), u_fn=u_fn,
                      t0=T0, record="full")

    # 동조 판정: 1:1 위상 잠금.
    #   구동 대비 진동자 peak 의 상대위상이 관측창 내내 표류하지 않아야 한다.
    #   주기 일치만 보면 세기 0 에서도 우연히 통과한다(고유주기 ~ 구동주기 인 열).
    ent = np.zeros(N, bool); obs = np.full(N, np.nan)
    drift = np.full(N, np.nan); alive = np.zeros(N, bool)
    rot = np.full(N, np.nan)
    span = 60 * tau
    for j in range(N):
        b = tr[:, PHASE_VAR, j]
        alive[j] = (b.max() - b.min()) / max(b.mean(), 1e-9) > 0.10   # 진폭 사멸 배제
        pk = peaks(b, t0=T0)
        if len(pk) < 10:
            continue
        obs[j] = np.median(np.diff(pk))
        # 회전수 rho = 진동자 주기 수 / 구동 주기 수. 1:1 동조면 rho ~= 1
        rot[j] = (len(pk) - 1) * T[j] / (pk[-1] - pk[0])
        phi = (pk % T[j]) / T[j]                    # 구동 대비 상대위상 [0,1)
        u = np.unwrap(phi * 2 * np.pi) / (2 * np.pi)
        cyc = (pk - pk[0]) / T[j]
        fit = np.polyfit(cyc, u, 1)
        drift[j] = abs(fit[0])                      # 주기당 위상 표류
        resid = u - np.polyval(fit, cyc)
        ent[j] = (alive[j] and abs(rot[j] - 1.0) < 0.02
                  and drift[j] < 0.004 and resid.std() < 0.02)

    return pd.DataFrame({"strength": S, "drive_period_h": T,
                         "period_ratio": T / tau, "observed_period_h": obs,
                         "rotation_number": rot, "phase_drift_per_cycle": drift,
                         "oscillating": alive, "entrained": ent})


def main():
    print("=" * 78)
    per, amp, rel, pm, _ = characterize(np.array([KDEG0]))
    tau = float(per[0])
    print(f"[1] 기저 진동:  주기 {tau:.3f} h   진폭 {amp[0]:.4f}   상대진폭 {rel[0]:.3f}")

    kdeg_sweep(tau)

    print("\n[3] 위상반응곡선 (PRC) — 중력 펄스 1 h")
    frames = []
    for label, kd, stg in [("1G baseline (약한 자극 0.05)", KDEG0, 0.05),
                           ("1G baseline (중간 자극 0.10)", KDEG0, 0.10),
                           ("microgravity-like kdeg80% (0.05)", KDEG0 * 0.8, 0.05)]:
        prc = compute_prc(tau, kdeg=kd, strength=stg)
        if prc is None:
            print(f"    {label}: 한계순환 없음 — 생략"); continue
        prc["condition"] = label
        prc["kdeg"] = kd
        prc["strength"] = stg
        frames.append(prc)
        v = prc.phase_shift_h
        adv, dly = prc.loc[v.idxmax()], prc.loc[v.idxmin()]
        ptt = v.max() - v.min()
        print(f"    {label}:")
        print(f"       최대 전진 {adv.phase_shift_h:+.3f} h @ CT {adv.CT_h:5.1f}   "
              f"최대 지연 {dly.phase_shift_h:+.3f} h @ CT {dly.CT_h:5.1f}")
        print(f"       PRC 진폭(peak-to-trough) {ptt:.3f} h  "
              f"→ {'type 1' if ptt < tau/2 else 'type 0'}")
        dead = prc[v.abs() < 0.05 * ptt]
        if len(dead):
            print(f"       무반응 구간(dead zone) 약 CT "
                  f"{dead.CT_h.min():.0f}-{dead.CT_h.max():.0f}")

    if frames:
        allprc = pd.concat(frames)
        allprc.to_csv(os.path.join(OUT, "model_prc.csv"), index=False, encoding="utf-8-sig")
        a = frames[0].phase_shift_h
        if len(frames) >= 3:
            c = frames[2].phase_shift_h
            print(f"\n[5] 미세중력 상태에서 PRC 진폭: "
                  f"{a.max()-a.min():.3f} h -> {c.max()-c.min():.3f} h "
                  f"({((c.max()-c.min())/(a.max()-a.min())-1)*100:+.1f}%)")
            print("    → 웻랩에서 원심분리 세기를 정할 때 직접 쓰이는 수치")

        # 텍스트 PRC
        print(f"\n    PRC ({frames[0].condition.iloc[0]}):")
        p0 = frames[0].sort_values("CT_h")
        sc = max(abs(p0.phase_shift_h).max(), 1e-9)
        for _, r in p0.iloc[::2].iterrows():
            pos = int(round(r.phase_shift_h / sc * 20))
            bar = " " * (20 + min(pos, 0)) + ("#" * abs(pos) if pos else "|")
            print(f"      CT {r.CT_h:5.1f} {r.phase_shift_h:+6.3f} |{bar}")

    print("\n[4] Arnold tongue — 펄스 세기 x 구동주기")
    at = arnold_tongue(tau)
    at.to_csv(os.path.join(OUT, "model_arnold_tongue.csv"), index=False, encoding="utf-8-sig")
    piv = at.pivot_table(index="strength", columns="period_ratio", values="entrained")
    cols = piv.columns.to_numpy()
    print("        " + "".join("v" if abs(c - 1.0) < 0.007 else " " for c in cols))
    for s, row in piv.iterrows():
        print(f"    {s:.3f} |" + "".join("#" if v else "." for v in row.to_numpy()) + "|")
    print(f"    열: 구동주기/고유주기 {cols.min():.2f} ~ {cols.max():.2f}")
    print("\n    세기별 1:1 동조 가능 주기 범위 (= 동조 대역폭):")
    for s in sorted(at.strength.unique())[::2]:
        e = at[(at.strength == s) & at.entrained]
        if len(e):
            print(f"      세기 {s:.3f}: {e.drive_period_h.min():.2f} ~ {e.drive_period_h.max():.2f} h "
                  f"(폭 {e.drive_period_h.max()-e.drive_period_h.min():.2f} h)")
        else:
            print(f"      세기 {s:.3f}: 동조 없음")
    print("\n    ※ 세기 0.20 초과 영역은 고차(n:m) 잠금과 진폭 사멸이 섞여")
    print("       1:1 동조 판정이 의미를 잃는다. 별도 회전수 분석이 필요하며 본 계획 범위 밖.")
    print(f"\n-> {OUT}/model_kdeg_sweep.csv, model_prc.csv, model_arnold_tongue.csv")


if __name__ == "__main__":
    main()
