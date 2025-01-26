FC = gfortran
TARGET = main
SRC = main.f90
FFLAGS := -O0 -fbacktrace -fcheck=all

#モジュールファイルの出力場所
MODDIR = ./mod

#----------------------------------------------------------
# 以下共通
#----------------------------------------------------------
# モジュールファイルの場所.
INCLUDE = ../../build/mod
# コンパイル済みアーカイブファイルの場所.
LDFLAGS = -L../../build/lib -lpartras

${TARGET}: ${SRC}
	${FC} ${SRC} -J${MODDIR} -I${INCLUDE} ${LDFLAGS} -o $@ ${FFLAGS}
