<#
.SYNOPSIS
  git push 가 자격증명 단계에서 실패할 때 재부팅 대신 시도하는 복구 스크립트.

.DESCRIPTION
  증상
    git add / commit / log 는 정상인데 push 만 아래 오류로 실패한다.

      *** fatal error - add_item ("\??\C:\Program Files\Git", "/", ...) failed, errno 1
      fatal: could not read Username for 'https://github.com': terminal prompts disabled

  원인 (확인된 부분)
    git 은 자격증명 헬퍼를 명령 문자열 하나로 만들어 셸에 넘긴다. GIT_TRACE 로 확인한 실제 호출:

      run_command:   'git credential-manager get'
      start_command: 'C:/Program Files/Git/usr/bin/sh.exe' -c 'git credential-manager get'

    문자열에 공백이 있으므로 항상 sh.exe 를 거친다. 어떤 헬퍼를 쓰든 '<헬퍼> get' 형태라
    공백을 피할 수 없다. 따라서 MSYS2 런타임이 불안정해지면 HTTPS push 만 골라서 실패한다.

  원인 (미확정)
    MSYS2 가 왜 간헐적으로 죽는지는 규명하지 못했다. 아래는 배제했다.
      - msys-2.0.dll 중복 설치 : 없음 (1개)
      - Mandatory ASLR         : 꺼져 있음
      - fstab 마운트 과다      : 항목 2개뿐
      - 환경변수 블록 과대     : 3.5KB, 정상
    남은 가설은 좀비 MSYS2 프로세스가 공유 프로세스 테이블을 오염시킨다는 것이다.
    재부팅으로 풀린다는 점과는 맞지만, 고장을 재현하지 못해 확정하지 못했다.
    이 스크립트는 그 가설에 기반한 시도이며, 듣지 않으면 재부팅이 여전히 필요하다.

.NOTES
  근본 해결은 SSH 로 전환하는 것이다. push 경로에서 MSYS2 를 완전히 제거한다.
  자세한 내용은 이 파일 맨 아래 주석 참조.
#>

[CmdletBinding()]
param(
    [switch]$Force,          # 확인 없이 좀비 프로세스 종료
    [switch]$Push            # 복구 후 git push origin main 까지 실행
)

$ErrorActionPreference = 'Continue'
$GitRoot = 'C:\Program Files\Git'
$Sh      = Join-Path $GitRoot 'usr\bin\sh.exe'

function Test-Shell {
    try {
        $out = & $Sh -c "echo ok" 2>&1
        return ($LASTEXITCODE -eq 0 -and "$out" -match 'ok')
    } catch { return $false }
}

Write-Host "`n[1/4] sh.exe 상태 확인" -ForegroundColor Cyan
if (Test-Shell) {
    Write-Host "  정상. MSYS2 런타임에는 문제가 없다." -ForegroundColor Green
    Write-Host "  push 가 실패한다면 원인이 다른 데 있다 (네트워크, 토큰 만료, 권한)."
} else {
    Write-Host "  실패. MSYS2 런타임이 깨져 있다." -ForegroundColor Yellow
}

Write-Host "`n[2/4] 좀비 MSYS2 프로세스 탐색" -ForegroundColor Cyan
# 핸들 0 개는 정상 셸에서 나올 수 없는 값이다. 죽다 만 프로세스로 본다.
$zombies = Get-Process -Name bash, sh, ssh-agent, mintty -ErrorAction SilentlyContinue |
           Where-Object { $_.HandleCount -eq 0 }
$msysAll = Get-Process -Name bash, sh, mintty -ErrorAction SilentlyContinue

Write-Host ("  MSYS2 프로세스 {0}개, 그중 좀비 {1}개" -f @($msysAll).Count, @($zombies).Count)
foreach ($z in $zombies) {
    Write-Host ("    PID {0}  {1}  핸들={2}  스레드={3}" -f $z.Id, $z.ProcessName, $z.HandleCount, $z.Threads.Count)
}

if (@($zombies).Count -eq 0) {
    Write-Host "  종료할 대상이 없다." -ForegroundColor Green
} else {
    $go = $Force
    if (-not $go) {
        $ans = Read-Host "  이 프로세스들을 종료하겠는가? (y/N)"
        $go = ($ans -eq 'y')
    }
    if ($go) {
        Write-Host "`n[3/4] 종료" -ForegroundColor Cyan
        foreach ($z in $zombies) {
            try { Stop-Process -Id $z.Id -Force -ErrorAction Stop; Write-Host ("  종료됨 PID {0}" -f $z.Id) -ForegroundColor Green }
            catch { Write-Host ("  종료 실패 PID {0}: {1}" -f $z.Id, $_.Exception.Message) -ForegroundColor Red }
        }
        Start-Sleep -Milliseconds 800
    } else {
        Write-Host "  건너뜀."
    }
}

Write-Host "`n[4/4] 복구 확인" -ForegroundColor Cyan
if (Test-Shell) {
    Write-Host "  sh.exe 정상 작동." -ForegroundColor Green
    if ($Push) {
        Write-Host "`n  git push origin main 실행" -ForegroundColor Cyan
        git push origin main
        if ($LASTEXITCODE -eq 0) { Write-Host "  push 성공." -ForegroundColor Green }
        else { Write-Host "  push 실패. 아래 '근본 해결' 참조." -ForegroundColor Red }
    } else {
        Write-Host "  이제 git push 를 시도하라. (-Push 를 주면 여기서 바로 실행한다)"
    }
} else {
    Write-Host "  여전히 실패." -ForegroundColor Red
    Write-Host @"

  이 저장소(KARA)는 2026-08-02 에 SSH 로 전환했으므로 이 오류의 영향을 받지 않는다.
  아래는 아직 HTTPS 를 쓰는 다른 저장소를 위한 안내다.

  (1) 재부팅 — 확실하지만 매번 해야 한다.

  (2) SSH 로 전환 — push 경로에서 MSYS2 를 완전히 제거한다. 권장.

      확인된 사실: 공백 없는 단일 경로를 core.sshCommand 로 주면
      git 이 sh.exe 를 거치지 않고 ssh.exe 를 직접 실행한다.
      KARA 에서 dry-run 으로 검증한 결과 sh.exe 호출 0 회였다.

        trace: start_command: C:/Windows/System32/OpenSSH/ssh.exe git@github.com 'git-receive-pack ...'
        (sh -c 래핑 없음)

      필요한 작업
        a. 공개키를 GitHub 에 등록 (https://github.com/settings/keys)
           키 파일: %USERPROFILE%\.ssh\id_ed25519.pub
           (이 키는 등록 완료. 암호가 없어 무인 push 가 가능하다)

        b. 설정 (a 를 마친 뒤에 실행할 것. 먼저 하면 push 가 아예 막힌다)
           git config --global core.sshCommand "C:/Windows/System32/OpenSSH/ssh.exe"   # 전역, 이미 적용됨
           git remote set-url origin git@github.com:<계정>/<저장소>.git

        c. 검증
           C:\Windows\System32\OpenSSH\ssh.exe -T git@github.com
           -> "Hi josh-min99! You've successfully authenticated" 가 나오면 성공

"@ -ForegroundColor Yellow
}

Write-Host ""
