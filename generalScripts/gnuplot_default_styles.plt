
SOLID=1
DASHEDBIG=2
DASHEDSMALL=3
DOTTED=4 

PLUS=1
EX=2
STAR=3
BOXEMPTY=4
BOXFILLED=5
CIRCLEEMPTY=6
CIRCLEFULL=7

RED=1
GREEN=2
BLACK=-1
BLUE=3
MAGENTA=4

#whatever styles 1, 2, 3 ... are will be used by default for the first, second, third etc data series if
#set style increment user is set

#red solid
    REDSOLID=1
	set style line REDSOLID lw 2  lt SOLID     lc RED pt PLUS       ps 1
#red dashed
    REDDASHED=5
	set style line REDDASHED lw 2 lt DASHEDBIG lc RED pt CIRCLEFULL ps 1
#red dotted
    REDDOTTED=9
	set style line REDDOTTED lw 3 lt DOTTED    lc RED pt EX         ps 1

#blue solid
    BLUESOLID=4
	set style line BLUESOLID lw 2  lt SOLID     lc BLUE pt PLUS       ps 1
#blue dashed
    BLUEDASHED=2
	set style line BLUEDASHED lw 2 lt DASHEDBIG lc BLUE pt CIRCLEFULL ps 1
#blue dotted
    BLUEDOTTED=6
	set style line BLUEDOTTED lw 3 lt DOTTED    lc BLUE pt EX         ps 1

#black solid
    BLACKSOLID=7
	set style line BLACKSOLID lw 2  lt SOLID     lc BLACK pt PLUS       ps 1
#black dashed
    BLACKDASHED=8
	set style line BLACKDASHED lw 2 lt DASHEDBIG lc BLACK pt CIRCLEFULL ps 1
#black dotted
    BLACKDOTTED=3
	set style line BLACKDOTTED lw 3 lt DOTTED    lc BLACK pt EX         ps 1


#large magenta circle
	#set style line 5 lw 2 ps 1.5 pt 65 lc 4
	#set style line 5 lw 2 ps 1.5 pt 6 lc 4
#green dashed
    #GREENDASHED=6
	#set style line 6 lw 2 lt 2 lc 2
#black ?
	#set style line 7 lw 2 ps 1.5 pt 6 lc -1
#green ?
	#set style line 8 lw 2 ps 1.5 pt 6 lc 2



