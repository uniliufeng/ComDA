! date.F90 - GSM 風日時ルーチン
! vi: set sw=2 ts=72:
!
! date モジュールは長さ 5 の整数型配列として保存した日時演算を提供する。
!
! ここで日時は以下のように格納される:
!
!       id(1): 西暦年
!       id(2): 月
!       id(3): 日
!       id(4): 時
!       id(5): 分*60+秒
!
! 特に、西暦年と月が 0 になっていることがある。このときの日の欄は
! 0 年 0 月 1 日 (つまり紀元前2年12月1日) を 1 とする通日である。
! 差を取るなどの演算では通日を用いるが、通日と普通の表記を区別しないで
! よいように上記紀元が選択されている。
!
! 注意: 通日の計算など、整数型が 16bit ではできない処理がある。

module date

! 暦法
!       暦法とは通日から年月を得る（あるいは逆演算の）方式である。
! 
!       360: 毎月 30 日の暦
!       365: 毎年 365 日の暦
!       36525: ユリウス暦
!       3652425: グレゴリオ・ユリウス暦自動切換 (1752 年 9 月 3 日)
!       その他の値: 3652425 とみなされる
!
  integer, save:: calendar_type = 3652425

contains

!=======================================================================
! date_diff - 日時の差
!       idate から idate_base を引いた差を idate_diff に格納する。
!       差は日、時、秒だけで表現され、時は -23 .. 23, 秒は -3599 .. 3599
!       date_calendar による暦法設定によって動作が変わる.

subroutine date_diff(idate, idate_base, idate_diff)
  integer, intent(in):: idate(5), idate_base(5)
  integer, intent(out):: idate_diff(5)
  integer:: now(5), base(5)
  now(:) = idate(:)
  base(:) = idate_base(:)
  call date_compact(now)
  call date_compact(base)
  idate_diff(:) = now(:) - base(:)
  call date_diff_normalize(idate_diff)
end subroutine

!=======================================================================
! date_diff_normalize - 日時の差をわかりやすい形式に書き換え
!       日付の差 id を表示に適した形式に書き換える。

subroutine date_diff_normalize(id)
  integer, intent(inout):: id(5)
  ! 年月は日時と可約でないのでこれでおしまい
  id(2) = id(1) * 12 + id(2) - 1
  id(1) = id(2) / 12
  id(2) = mod(id(2), 12) + 1
  ! 日時を日と秒に分離する (全部秒にすると数十年で溢れるため)
  id(3) = id(3) + id(4) / 24 + id(5) / 86400
  id(5) = mod(id(5), 86400) + mod(id(4), 24)
  id(4) = 0
  ! もし日と秒の符号が逆なら
  if (id(5) < 0 .and. id(3) > 0) then
    id(3) = id(3) - 1
    id(5) = id(5) + 86400
  else if (id(5) > 0 .and. id(3) < 0) then
    id(3) = id(3) + 1
    id(5) = id(5) - 86400
  endif
  ! 秒から時を生成
  id(4) = id(5) / 3600
  id(5) = mod(id(5), 3600)
end subroutine

!=======================================================================
! date_compare - 日時の前後比較
!       日付 id, id_base を比較し、判定結果を cmp に書き込む。
!       cmp = 1: date が後, cmp = 0: 同時, cmp = -1: date が前
!       周期性は意味がないので無視。

subroutine date_compare(id, id_base, cmp)
  integer, intent(in):: id(5), id_base(5)
  integer, intent(out):: cmp
  integer:: id_diff(5), i
  call date_diff(id, id_base, id_diff)
  do, i = 1, 5
    cmp = sign(id_diff(i), id_diff(i))
    if (cmp /= 0) exit
  enddo
end subroutine

!=======================================================================
! date_diff_compact - 日時の差を通月・通日・通秒に書き換え
!       日付の差 id を演算に適した形式に書き換える。

subroutine date_diff_compact(id)
  integer, intent(inout):: id(5)
  id(3) = id(3) + id(4) / 24 + id(5) / 86400 
  id(5) = mod(id(5), 86400) + mod(id(4), 24) * 3600
  id(4) = 0
  id(2) = id(2) + id(1) * 12 
  id(1) = 0
end subroutine

!=======================================================================
! date_diff_div - 日時の差を除算
!       日時間隔 a を日時間隔 b で割った結果 f を求める。
!       日時間隔は年月と日時の単位を混合してはならない。

subroutine date_diff_div(a, b, f)
  integer, intent(in):: a(5), b(5)
  double precision:: f
  integer:: ia(5), ib(5)
  ia(:) = a(:)
  call date_diff_compact(ia)
  ib(:) = b(:)
  call date_diff_compact(ib)
  if (ib(2) /= 0) then
    f = dble(ia(2)) / dble(ib(2))
  else if (ib(3) /= 0) then
    f = dble(ia(3)) + dble(ia(5)) / 86400.0
    f = f / (dble(ib(3)) + dble(ib(5)) / 86400.0)
  else if (ib(5) /= 0) then
    f = (ia(5) + ia(3) * 86400.0d0) / dble(ib(5))
  else
    print *, "date_diff_div: division by zero, a=", a
    f = 0.0d0
  endif
end subroutine

!=======================================================================
! date_modulo - 時間軸上の日時をある周期で折り返す
!       日時 id が id_origin から id_cycle の範囲に収まるように
!       id_cycle の整数倍を加減する。

subroutine date_modulo(id, id_origin, id_cycle)
  integer, intent(inout):: id(5)
  integer, intent(in):: id_origin(5), id_cycle(5)
  integer:: idiff(5), idiff_cycle(5)
  double precision:: ratio
  call date_diff(id, id_origin, idiff)
  call date_diff(id_origin + id_cycle, id_origin, idiff_cycle)
  call date_diff_div(idiff, idiff_cycle, ratio)
  if (ratio >= 0.0d0 .and. ratio < 1.0d0) goto 1000
  id = id - id_cycle * floor(ratio)
  1000 continue
  call date_normalize(id)
end subroutine

!=======================================================================
! date_calendar - 暦法を設定する
!       caltype == 0 で呼ぶと問い合わせ
!       それ以外の場合は設定。意味は calendar_type のコメント参照

subroutine date_calendar(caltype)
  integer, intent(inout):: caltype
  if (caltype == 0) then
    caltype = calendar_type
  else
    calendar_type = caltype
  endif
end subroutine

!=======================================================================
! date_compact - 日時表記を通日・秒に変換
!       idate を紀元前1年3月あたりを0とする通日、日内通秒 だけで表現し直す。
!       date_calendar による暦法設定によって動作が変わる.

subroutine date_compact(idate)
  integer, intent(inout):: idate(5)
  integer, parameter:: four_years = 365 * 4 + 1
  integer:: second, hour, iday_sec, iday_hour, year, month, century
  ! 時・秒を日と通秒に分離。通秒は負にしない
  second = modulo(idate(5), 86400)
  hour = modulo(idate(4), 24)
  iday_sec = (idate(5) - second) / 86400
  iday_hour = (idate(4) - hour) / 24
  idate(5) = second + hour * 3600
  idate(3) = idate(3) + iday_sec + iday_hour
  idate(4) = 0
  ! 1月と2月は前年の扱いにする
  month = modulo(idate(2) - 3, 12) + 3
  year = idate(1) + (idate(2) - month) / 12
  century = (year - modulo(year, 100)) / 100 + 1
  idate(1:2) = 0
  if (calendar_type == 360) then
    ! 1 年が 360 日の暦
    idate(3) = idate(3) + (year * 12 + month) * 30
    return
  endif
  ! 月日を毎年3月1日を 1 とする通日に変換
  idate(3) = idate(3) + (month * 306 - 914) / 10
  ! 年と年内通日から世紀通日を算出
  if (calendar_type == 365) then
    ! 1 年が 365 日の暦
    idate(3) = idate(3) + year * 365 + 90
    return
  else
    idate(3) = idate(3) + (year * four_years - modulo(year * four_years, 4))/ 4
    ! デフォルトでのグレゴリオ暦移行の日付は英国のものである
    if (calendar_type == 36525 .or. idate(3) < 640116) then
      ! ユリウス暦
      idate(3) = idate(3) + 91
    else
      ! 世紀による補正: グレゴリオ暦, 通日はユリウス暦と互換
      idate(3) = idate(3) - (century * 3 - modulo(century * 3, 4)) / 4 + 93
    endif
  endif
end subroutine

!=======================================================================
! date_normalize - 通日・通秒や加算結果などによる不整日時を直す
!       date_calendar による暦法設定によって動作が変わる.

subroutine date_normalize(idate)
  integer, intent(inout):: idate(5)
  integer:: second, hour, iday_sec, iday_hour, day, year, month
  integer, parameter:: four_year = 365 * 4 + 1
  integer, parameter:: four_century = 365 * 400 + 97

  ! 入力が通日になっていなければ通日変換
  if (idate(1) /= 0 .or. idate(2) /= 0) then
    call date_compact(idate)
  endif

  ! 通秒から時秒を分離。通秒は負にしない
  second = modulo(idate(5), 86400)
  hour = modulo(idate(4), 24)
  iday_sec = (idate(5) - second) / 86400
  iday_hour = (idate(4) - hour) / 24
  idate(5) = modulo(second, 3600)
  idate(4) = hour + second / 3600
  idate(3) = idate(3) + iday_sec + iday_hour

  ! 暦法によって日から年月の算出は異なる
  if (calendar_type == 360) then
    day = modulo(idate(3) - 31, 360)
    idate(1) = (idate(3) - 31 - day) / 360
    idate(2) = day / 30 + 1
    idate(3) = modulo(day, 30) + 1
    return
  endif
  if (calendar_type == 365) then
    day = modulo(idate(3) - 91, 365)
    idate(1) = (idate(3) - 91 - day) / 365
  else
    ! オリンピックのある年の3月1日を0とする通日を計算
    if (calendar_type == 36525 .or. idate(3) < 640196) then
      ! ユリウス暦
      day = modulo(idate(3) - 92, four_year)
      year = (idate(3) - 92 - day) / four_year * 4
    else
      ! グレゴリオ暦
      day = modulo(idate(3) - 94, four_century)
      year = (idate(3) - 94 - day) / four_century * 400
      if (day == four_century - 1) then
        year = year + 300
        day = 36524
      else
        year = year + day / 36524 * 100
        day = modulo(day, 36524)
      endif
      year = year + day / four_year * 4
      day = modulo(day, four_year)
    endif
    ! オリンピックのある年の3月1日を0とする通日から
    ! 3月ではじまる年と年内通日を計算
    if (day == four_year - 1) then
      idate(1) = year + 3
      day = 365
    else
      idate(1) = year + day / 365
      day = modulo(day, 365)
    endif
  endif
  ! 3月1日を0とする通日から月日を計算し、年を割礼年始に修正
  day = day * 10 + 922
  month = day / 306
  idate(2) = mod(month - 1, 12) + 1
  idate(1) = idate(1) + (month - idate(2)) / 12
  idate(3) = mod(day, 306) / 10 + 1
end subroutine

!=======================================================================
end module
