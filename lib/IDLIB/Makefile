include ../make.inc
#
# Specify the Fortran compiler FORTRAN and its flags FFLAGS.
#

#FC=lf95
#FFLAGS=--o2

#FC=fort77
#FFLAGS=-O3 -!P

#FC=gfortran
FFLAGS=-O3 -w

#
# Specify the locations BLASLIB and LAPACK_LIB
# of the BLAS and LAPACK archives (libraries).
#

#BLASLIB=/home/tygert/matapprox/id/lbin/blas_LINUX.a
#LAPACK_LIB=/home/tygert/matapprox/id/lbin/lapack_LINUX.a

#BLASLIB=~/matapprox/id/bin/blas_LINUX.a
#LAPACK_LIB=~/matapprox/id/bin/lapack_LINUX.a

BLASLIB=../librefblas.a
LAPACK_LIB=../liblapack.a

#
# Specify the archiver ARCH and its flags ARCHFLAGS used
# to create the archive (library) .a file.
#

#ARCH=ar
#ARCHFLAGS=cr

#
# Set RANLIB=echo if your system does not have ranlib.
#

#RANLIB=ranlib

###############################################################################

ALL_ID_OBJS=prini.o dfft.o id_rand.o id_rtrans.o \
            idd_snorm.o idd_sfft.o idd_frm.o \
            idd_house.o idd_qrpiv.o idd_id.o idd_svd.o idd_id2svd.o \
            iddr_rid.o iddr_rsvd.o iddp_rid.o iddp_rsvd.o \
            iddr_aid.o iddr_asvd.o iddp_aid.o iddp_asvd.o \
            idz_snorm.o idz_sfft.o idz_frm.o \
            idz_house.o idz_qrpiv.o idz_id.o idz_svd.o idz_id2svd.o \
            idzr_rid.o idzr_rsvd.o idzp_rid.o idzp_rsvd.o \
            idzr_aid.o idzr_asvd.o idzp_aid.o idzp_asvd.o

ID_RAND=id_rand.o prini.o
ID_RTRANS=id_rtrans.o prini.o id_rand.o
IDD_SNORM=idd_snorm.o prini.o id_rand.o
IDD_SFFT=idd_sfft.o prini.o dfft.o id_rand.o
IDD_HOUSE=idd_house.o prini.o
IDD_QRPIV=idd_qrpiv.o prini.o idd_house.o
IDD_ID=idd_id.o prini.o idd_house.o idd_qrpiv.o
IDD_SVD=idd_svd.o prini.o idd_house.o idd_qrpiv.o id_rand.o $(LAPACK_LIB) \
        $(BLASLIB)
IDD_ID2SVD=idd_id2svd.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
           $(LAPACK_LIB) $(BLASLIB)
IDD_FRM=idd_frm.o prini.o dfft.o id_rand.o id_rtrans.o idd_house.o \
        idd_qrpiv.o idd_id.o idd_sfft.o idd_svd.o $(LAPACK_LIB) $(BLASLIB)
IDDR_RID=iddr_rid.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o
IDDR_RSVD=iddr_rsvd.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
          iddr_rid.o idd_id2svd.o $(LAPACK_LIB) $(BLASLIB)
IDDP_RID=iddp_rid.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o
IDDP_RSVD=iddp_rsvd.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
          iddp_rid.o idd_id2svd.o $(LAPACK_LIB) $(BLASLIB)
IDDR_AID=iddr_aid.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
         idd_sfft.o id_rtrans.o idd_frm.o dfft.o
IDDR_ASVD=iddr_asvd.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
          idd_sfft.o id_rtrans.o idd_frm.o iddr_aid.o idd_id2svd.o \
          dfft.o $(LAPACK_LIB) $(BLASLIB)
IDDP_AID=iddp_aid.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
         idd_sfft.o id_rtrans.o idd_frm.o dfft.o
IDDP_ASVD=iddp_asvd.o prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o \
          idd_sfft.o id_rtrans.o idd_frm.o iddp_aid.o idd_id2svd.o \
          dfft.o $(LAPACK_LIB) $(BLASLIB)
IDZ_SNORM=idz_snorm.o prini.o id_rand.o
IDZ_SFFT=idz_sfft.o prini.o dfft.o id_rand.o
IDZ_HOUSE=idz_house.o prini.o
IDZ_QRPIV=idz_qrpiv.o prini.o idz_house.o
IDZ_ID=idz_id.o prini.o idz_house.o idz_qrpiv.o
IDZ_SVD=idz_svd.o prini.o idz_house.o idz_qrpiv.o id_rand.o $(LAPACK_LIB) \
        $(BLASLIB)
IDZ_ID2SVD=idz_id2svd.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
           $(LAPACK_LIB) $(BLASLIB)
IDZ_FRM=idz_frm.o prini.o dfft.o id_rand.o id_rtrans.o idz_house.o \
        idz_qrpiv.o idz_id.o idz_sfft.o idz_svd.o $(LAPACK_LIB) $(BLASLIB)
IDZR_RID=idzr_rid.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o
IDZR_RSVD=idzr_rsvd.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
          idzr_rid.o idz_id2svd.o $(LAPACK_LIB) $(BLASLIB)
IDZP_RID=idzp_rid.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o
IDZP_RSVD=idzp_rsvd.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
          idzp_rid.o idz_id2svd.o $(LAPACK_LIB) $(BLASLIB)
IDZR_AID=idzr_aid.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
         idz_sfft.o id_rtrans.o idz_frm.o dfft.o
IDZR_ASVD=idzr_asvd.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
          idz_sfft.o id_rtrans.o idz_frm.o idzr_aid.o idz_id2svd.o \
          dfft.o $(LAPACK_LIB) $(BLASLIB)
IDZP_AID=idzp_aid.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
         idz_sfft.o id_rtrans.o idz_frm.o dfft.o
IDZP_ASVD=idzp_asvd.o prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o \
          idz_sfft.o id_rtrans.o idz_frm.o idzp_aid.o idz_id2svd.o \
          dfft.o $(LAPACK_LIB) $(BLASLIB)

IDD_R=prini.o idd_house.o idd_qrpiv.o idd_id.o id_rand.o iddr_rid.o \
      iddp_rid.o idd_id2svd.o iddr_rsvd.o iddp_rsvd.o idd_snorm.o \
      $(LAPACK_LIB) $(BLASLIB)
IDZ_R=prini.o idz_house.o idz_qrpiv.o idz_id.o id_rand.o idzr_rid.o \
      idzp_rid.o idz_id2svd.o idzr_rsvd.o idzp_rsvd.o idz_snorm.o \
      $(LAPACK_LIB) $(BLASLIB)
IDD_A=prini.o idd_house.o idd_qrpiv.o idd_id.o idd_svd.o id_rand.o idd_sfft.o \
      id_rtrans.o idd_frm.o iddr_aid.o iddp_aid.o idd_id2svd.o iddr_asvd.o \
      iddp_asvd.o idd_snorm.o dfft.o $(LAPACK_LIB) $(BLASLIB)
IDZ_A=prini.o idz_house.o idz_qrpiv.o idz_id.o idz_svd.o id_rand.o idz_sfft.o \
      id_rtrans.o idz_frm.o idzr_aid.o idzp_aid.o idz_id2svd.o idzr_asvd.o \
      idzp_asvd.o idz_snorm.o dfft.o $(LAPACK_LIB) $(BLASLIB)

all: $(ALL_ID_OBJS) ../id_lib.a

clean:
	rm -f *.o

dfft:
	dfft is not supposed to be a target for make.
prini:
	prini is not supposed to a target for make.
id_rand:
	id_rand should be id_rand_test.
id_rtrans:
	id_rtrans should be id_rtrans_test.
idd_snorm:
	idd_snorm should be idd_snorm_test.
idd_sfft:
	idd_sfft should be idd_sfft_test.
idd_house:
	idd_house should be idd_house_test.
idd_qrpiv:
	idd_qrpiv should be idd_qrpiv_test.
idd_id:
	idd_id should be idd_id_test.
idd_svd:
	idd_svd should be idd_svd_test.
idd_id2svd:
	idd_id2svd should be idd_id2svd_test.
idd_frm:
	idd_frm should be idd_frm_test.
iddr_rid:
	iddr_rid should be iddr_rid_test.
iddr_rsvd:
	iddr_rsvd should be iddr_rsvd_test.
iddp_rid:
	iddp_rid should be iddp_rid_test.
iddp_rsvd:
	iddp_rsvd should be iddp_rsvd_test.
iddr_aid:
	iddr_aid should be iddr_aid_test.
iddr_asvd:
	iddr_asvd should be iddr_asvd_test.
iddp_aid:
	iddp_aid should be iddp_aid_test.
iddp_asvd:
	iddp_asvd should be iddp_asvd_test.
idz_snorm:
	idz_snorm should be idz_snorm_test.
idz_sfft:
	idz_sfft should be idz_sfft_test.
idz_house:
	idz_house should be idz_house_test.
idz_qrpiv:
	idz_qrpiv should be idz_qrpiv_test.
idz_id:
	idz_id should be idz_id_test.
idz_svd:
	idz_svd should be idz_svd_test.
idz_id2svd:
	idz_id2svd should be idz_id2svd_test.
idz_frm:
	idz_frm should be idz_frm_test.
idzr_rid:
	idzr_rid should be idzr_rid_test.
idzr_rsvd:
	idzr_rsvd should be idzr_rsvd_test.
idzp_rid:
	idzp_rid should be idzp_rid_test.
idzp_rsvd:
	idzp_rsvd should be idzp_rsvd_test.
idzr_aid:
	idzr_aid should be idzr_aid_test.
idzr_asvd:
	idzr_asvd should be idzr_asvd_test.
idzp_aid:
	idzp_aid should be idzp_aid_test.
idzp_asvd:
	idzp_asvd should be idzp_asvd_test.

../id_lib.a:
	$(ARCH) $(ARCHFLAGS) $@ $(ALL_ID_OBJS)
	$(RANLIB) $@

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< -o $@
