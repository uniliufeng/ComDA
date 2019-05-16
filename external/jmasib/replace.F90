! replace.F90 - 文字列関係の言語機能を補完するサブルーチン
! vi: set ts=70 sw=2:
! 2000-06-01 TOYODA Eizi

!=======================================================================
! uppercase - 小文字を大文字に置換する
! 文字型変数 string の中に英小文字があれば対応する英大文字に置換する。

subroutine uppercase(string)
  character(len = *), intent(inout):: string
  character(len = *), parameter:: UC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len = *), parameter:: LC = 'abcdefghijklmnopqrstuvwxyz'
  integer:: i, idx
  do, i = 1, len(string)
    idx = index(LC, string(i:i))
    if (idx == 0) cycle
    string(i:i) = UC(idx:idx)
  enddo
end subroutine

!=======================================================================
! getword - 文字列の先頭の語を抜き出す
! 文字型変数 string を空白で区切り、最初の語を word に格納し、
! 次の語以降を string に格納する。語がなければ word は空白となる

subroutine getword(string, word)
  character(len = *), intent(inout):: string
  character(len = *), intent(out):: word
  integer:: idx, strlen
  strlen = len(string)

  ! 先頭の空白の除去
  idx = verify(string, ' ')
  if (idx == 0) then
    word = ''
    return
  endif
  string = string(idx: strlen)

  ! 先頭の語の抜き出し
  idx = scan(string, ' ')
  if (idx == 0) then
    word = string
    string = ''
    return
  endif
  word = string(1: idx-1)
  string(1: idx-1) = ' '
  idx = verify(string, ' ')
  if (idx == 0) idx = 1
  string = string(idx:strlen)
end subroutine

!=======================================================================
! replace_int - 文字列中に整数値を埋め込む
! 文字型変数 string の中に pattern という部分文字列があれば、整数 int を
! 10 進表現したものに置き換える。pattern がみつからなければ何も起こらない。
! pattern の桁数が大きすぎる場合は右詰で左に 0 が入る。

subroutine replace_int(string, pattern, int)
  character(len = *), intent(inout):: string
  character(len = *), intent(in):: pattern
  integer, intent(in):: int
  character(len = 32):: format
  integer:: head, length, tail
  intrinsic index, trim, len_trim
continue
  head = index(string, trim(pattern))
  if (head < 1) return
  length = len_trim(pattern)
  tail = head + length - 1
  write(format, "('(i', i4, '.', i4, ')')") length, length
  write(string(head:tail), format) int
end subroutine

!=======================================================================
! replace_char - 文字列を別の文字列で置き換える
! 文字型変数 string の中に pattern という部分文字列があれば、char に
! 置き換える。pattern がみつからなければ何も起こらない。
! pattern と char の長さが違う場合は残りの部分が適切にシフトされる。

subroutine replace_char(string, pattern, char)
  character(len = *), intent(inout):: string
  character(len = *), intent(in):: pattern
  character(len = *), intent(in):: char
  integer:: head, shift, tail, i
  intrinsic index, trim, len_trim
continue
  head = index(string, trim(pattern))
  if (head < 1) return
  shift = len_trim(char) - len_trim(pattern)
  tail = head + len_trim(char) - 1
  if (shift > 0) then
    do i = len(string), tail, -1
      string(i:i) = string(i-shift:i-shift)
    enddo
  else if (shift < 0) then
    do i = tail, len(string) + shift, 1
      string(i:i) = string(i-shift:i-shift)
    enddo
    string(len(string) + shift: len(string)) = ''
  endif
  string(head: tail) = trim(char)
end subroutine

!/* for debugging */
#if DEBUG_MAIN
program main
  character(len = 60):: buffer
  character(len = 60):: homedir
  buffer = '$HOME/data/test%YYY_%M_%D_%H'
  call getenv('HOME', 4, homedir, len(homedir))
  print *, buffer
  call replace_char(buffer, '$HOME', homedir)
  print *, buffer
  call replace_char(buffer, 'data', 'dat')
  print *, buffer
  call replace_int(buffer, '%YYY', 2000)
  print *, buffer
  call replace_int(buffer, '%M', 5)
  print *, buffer
  call replace_int(buffer, '%D', 25)
  print *, buffer
  call replace_int(buffer, '%H', 23)
  print *, buffer
end program
#endif
