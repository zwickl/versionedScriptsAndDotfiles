" Vim 16-color scheme

hi clear
if exists("syntax_on")
    syntax reset
endif
let g:colors_name="caribbean"

if exists("g:caribbean_light") && g:caribbean_light
    set background=light
    hi Comment          ctermfg=Brown                                  cterm=NONE
    hi Identifier       ctermfg=DarkRed                                cterm=NONE
    hi Function         ctermfg=DarkRed                                cterm=NONE
    hi Statement        ctermfg=DarkBlue                               cterm=NONE
    hi Type             ctermfg=DarkBlue                               cterm=NONE
    hi Preproc          ctermfg=DarkBlue                               cterm=NONE
    hi Folded           ctermfg=DarkGray        ctermbg=Black          cterm=reverse,bold
    hi LineNr           ctermfg=DarkGray        ctermbg=LightGray      cterm=bold
    hi StatusLine       ctermfg=White           ctermbg=DarkBlue       cterm=bold
    hi StatusLineNC     ctermfg=DarkCyan        ctermbg=DarkBlue       cterm=NONE
    hi VertSplit        ctermfg=DarkBlue        ctermbg=DarkBlue       cterm=NONE
    hi TabLine          ctermfg=LightGray       ctermbg=DarkBlue       cterm=bold,underline
    hi TabLineFill      ctermfg=LightGray       ctermbg=DarkBlue       cterm=bold,underline
    hi TabLineSel       ctermfg=LightGray       ctermbg=Black          cterm=bold
else
    set background=dark
    hi Comment          ctermfg=DarkGray                               cterm=bold
    hi Identifier       ctermfg=Brown                                  cterm=NONE
    hi Function         ctermfg=Brown                                  cterm=NONE
    hi Statement        ctermfg=DarkCyan                               cterm=NONE
    hi Type             ctermfg=DarkCyan                               cterm=NONE
    hi Preproc          ctermfg=DarkCyan                               cterm=NONE
    hi Folded           ctermfg=DarkGray        ctermbg=Black          cterm=reverse,bold
    hi LineNr           ctermfg=DarkBlue
    hi StatusLine       ctermfg=White           ctermbg=DarkBlue       cterm=NONE
    hi StatusLineNC     ctermfg=LightBlue       ctermbg=DarkBlue       cterm=NONE
    hi VertSplit        ctermfg=DarkBlue        ctermbg=DarkBlue       cterm=NONE
    hi TabLine          ctermfg=DarkGray        ctermbg=DarkBlue       cterm=bold,underline
    hi TabLineFill      ctermfg=DarkGray        ctermbg=DarkBlue       cterm=bold,underline
    hi TabLineSel       ctermfg=Black           ctermbg=LightGray      cterm=NONE
endif

hi Constant             ctermfg=Magenta                                cterm=bold
hi String               ctermfg=DarkGreen                              cterm=NONE
hi Directory            ctermfg=DarkBlue                               cterm=NONE
hi Error                ctermfg=Red                                    cterm=bold,underline
hi ModeMsg              ctermfg=Black           ctermbg=Green          cterm=NONE
hi NonText              ctermfg=DarkBlue                               cterm=NONE
hi Search               ctermfg=Yellow          ctermbg=Black          cterm=bold,reverse
hi IncSearch            ctermfg=Yellow          ctermbg=Black          cterm=bold,reverse
hi Special              ctermfg=LightGreen                             cterm=NONE
" hi Visual               ctermfg=Black           ctermbg=Green          cterm=NONE
hi Visual               ctermfg=Green           ctermbg=Black          cterm=reverse
hi WarningMsg           ctermfg=DarkRed                                cterm=NONE

