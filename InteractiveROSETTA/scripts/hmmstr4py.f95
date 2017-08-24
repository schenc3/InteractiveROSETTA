module hmmstr
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::hmmstr_gamma
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::idatom
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::crossprod
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::read_profile
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::read_seq
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::readmodel
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::get_outtrans
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::get_intrans
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::flagrat_brdc_prod
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::bprof_bg
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::ramatype
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::getalpha
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::getbeta
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::read_backbone
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::release_model
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::set_flags
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !public::write_gca
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::get_angles
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::cal_tors
  
  private
   !------- module parameters -----------------------------------
   integer,parameter::maxseqlen=5000
   integer,parameter::m=20
   integer,parameter::mrama=11
   integer,parameter::mdssp=6
   integer,parameter::mctxt=10
   integer,parameter::smin=3
   integer,parameter::smax=27
   real,parameter::ncount=1.5
   real,parameter::eps=0.01
   real,parameter::alphaeps=0.01
   !------- module defined types --------------------------------
   !-------------------------------------------------------------
   ! HMMSTR model structure
   !-------------------------------------------------------------
   type model
     integer::N                          !! number of states
     real::bground(m)                    !! b background
     real::rground(mrama)                !! r background
     real::dground(mdssp)                !! d background
     real::cground(mctxt)                !! c background
     real,dimension(:),allocatable::prior
     real,dimension(:,:),allocatable::dih
     character(len=1),dimension(:),allocatable::ss
     real,dimension(:,:),allocatable::a   !! transitions
     real,dimension(:,:),allocatable::b   !! AA emissions
     real,dimension(:,:),allocatable::r   !! bb ang emissions
     real,dimension(:,:),allocatable::d   !! SS emissions
     real,dimension(:,:),allocatable::c   !! context emissions
     real,dimension(:,:),allocatable::logb
     character(len=1),dimension(6)::ss_map
   end type 
   !-------------------------------------------------------------
   ! Backbone coordinates
   !-------------------------------------------------------------
   type bb
     real,dimension(3,4)::xyz
   end type
   !-------------------------------------------------------------
   ! State-state transition probabilities.
   ! 0th value is th4e length of the array.
   !-------------------------------------------------------------
   type trans
    integer,dimension(:),allocatable::c
   end type
   !-------- module global variable -----------------------------
   type(model)::modelR
   type(trans),dimension(:),allocatable::intrans,outtrans
   real::sumd(mdssp),sumr(mrama),sumc(mctxt),gsum
   real,dimension(:,:),allocatable::calpha
   type(bb),dimension(:),allocatable::bbxyz
 contains
!------------------------------------------------------------------------------------------------------------
  !subroutine hmmstr_gamma(xgamma,nres,mfile,seqfile,pfile,use_str,pdbfile,chain,gca)
  subroutine hmmstr_gamma(nres,mfile,seq,backbone,gprofile)
   !---------------------------------------------
   ! Calculate the forward/backward algorithm
   ! and return the a postereiori probabilities
   ! of all HMM states at all positions.
   !---------------------------------------------
   implicit none
   !---------- dummy args -----------------------
   character(len=1000000),intent(in)::mfile
   integer,intent(inout)::nres
   real,dimension(:,:),intent(in)::backbone
   real,dimension(nres,20),intent(out)::gprofile
   integer,dimension(nres),intent(in)::seq
     !---------- subroutine variables -------------
   character(len=1),dimension(:),allocatable::ramaseq
   integer::t,nnode
   integer::i
   integer::iatom
   integer::ires
   integer,dimension(:),allocatable::iflag
   real,dimension(:),allocatable::ct
   real,dimension(:,:),allocatable::profile
   real,dimension(:,:),allocatable,target::angles
   real,dimension(:,:),allocatable::alpha
   real,dimension(:,:),allocatable::beta
   integer::inode
   real,dimension(282,nres)::xgamma
     !----------- read HMMSTR model ---------------
   call readmodel(mfile,sumd,sumr,sumc)
   call get_intrans
   call get_outtrans
   !----------------------------------------------
   !call readmodel(mfile,modelR,sumd,sumr,sumc)
   !call get_intrans(intrans,modelR)
   !call get_outtrans(outtrans,modelR)
   !------------ read sequence, structure data --
   if(allocated(bbxyz))    deallocate(bbxyz)
   if(allocated(profile))  deallocate(profile)
   !if(allocated(backbone)) deallocate(backbone)
   if(allocated(calpha))   deallocate(calpha)
   if(allocated(angles))   deallocate(angles)
   !if(allocated(resseq))   deallocate(resseq)
   !call read_seq(seqfile,seq,nres)
   !allocate(resseq(nres))
   allocate(calpha(3,nres));calpha=999.
   allocate(angles(3,nres));angles=0.
   !allocate(backbone(3,4*nres))
   allocate(bbxyz(nres))
   allocate(profile(20,nres))
   !if(use_str.and.pdbfile/=" ")then
   if(.True.)then
     !call read_backbone(pdbfile,chain,bbxyz,seq,nres,resseq)
     !call read_backbone(pdbfile,chain,seq,nres,resseq)
     !--------------------------------------------------------
     iatom=1
     do ires=1,nres
       bbxyz(ires)%xyz(:,1:4) = backbone(:,iatom:iatom+3)
       iatom = iatom+4
     enddo
     !--------------------------------------------------------
     do ires=1,nres
       calpha(1:3,ires)=bbxyz(ires)%xyz(1:3,2)
     enddo
     call get_angles(backbone,nres,angles)
   endif
   !if(pfile/=" ")then
   !  profile=0.
   !  call read_profile(pfile,profile,seq,nres)
   !else
   profile=0.
   do ires=1,nres
   	profile(seq(ires),ires)=1.00
   enddo
   !endif
   write(*,'(a,20f8.4)')"Background AA ",modelR%bground
   write(*,'(a,11f8.4)')"Background RAMA ",modelR%rground
   write(*,'(a,6f8.4)')"Background DSSP ",modelR%dground
   write(*,'(a,10f8.4)')"Background CTXT ",modelR%cground
   if(allocated(iflag))   deallocate(iflag)
   if(allocated(alpha))   deallocate(alpha)
   if(allocated(ct))      deallocate(ct)
   if(allocated(beta))    deallocate(beta)
   !if(allocated(xgamma))  deallocate(xgamma)
   if(allocated(ramaseq)) deallocate(ramaseq)
   allocate(ct(nres))
   allocate(iflag(nres))
   allocate(ramaseq(nres))
   allocate(alpha(modelR%n,nres))
   allocate(beta(modelR%n,nres))
   !allocate(xgamma(modelR%n,nres))
   !call set_flags(ramaseq,nres,use_str,angles,iflag)
   call set_flags(ramaseq,nres,angles,iflag)
   call getalpha(ct,alpha,nres,profile,ramaseq,iflag)
   call getbeta(ct,beta,nres,profile,ramaseq,iflag)
   !call getalpha(ct,alpha,nres,modelR,profile,intrans,ramaseq,iflag)
   !call getbeta(ct,beta,nres,modelR,profile,outtrans,ramaseq,iflag)
   nnode = modelR%N
   do t=1,nres
     do i=1,modelR%N
       xgamma(i,t)=alpha(i,t)*beta(i,t)/ct(t)
     enddo
     !!Write GPROFILE line
     gprofile(t,1:20)=0.
     do inode=1,modelR%N
       gprofile(t,1:20)=gprofile(t,1:20)+xgamma(inode,t)*modelR%b(inode,1:20)
     enddo
     gprofile(t,1:20)=gprofile(t,1:20)/sum(gprofile(t,1:20))
   enddo
   
   !if (present(gca)) &
     !call write_gca(gca,mfile,nres,seq,angles,calpha,modelR%N,xgamma,modelR,profile,resseq,ramaseq)
   !  call write_gca(gca,mfile,nres,seq,angles,calpha,nnode,xgamma,profile,resseq,ramaseq) 
   !!------------------------------------------------------------------------------------------------
   if(allocated(bbxyz))    deallocate(bbxyz)
   if(allocated(iflag))    deallocate(iflag)
   if(allocated(alpha))    deallocate(alpha)
   if(allocated(ct))       deallocate(ct)
   if(allocated(beta))     deallocate(beta)
   if(allocated(ramaseq))  deallocate(ramaseq)
   if(allocated(profile))  deallocate(profile)
   !if(allocated(backbone)) deallocate(backbone)
   if(allocated(angles))   deallocate(angles)
   !------------------------------------------------------------------------------------------------
   call release_model
   if(allocated(intrans))    deallocate(intrans)
   if(allocated(outtrans))   deallocate(outtrans)
  end subroutine hmmstr_gamma
!-------------------------------------------------------------------------------------------------------------
  !subroutine write_gca(gcafile,mfile,nres,seq,angles,calpha,nnode,xgamma,modelR,profile,resseq,ramaseq)
  !subroutine write_gca(gcafile,mfile,nres,seq,angles,calpha,nnode,xgamma,profile,resseq,ramaseq) 
  !  character(len=*),intent(in)::gcafile,mfile
  !  character(len=1),dimension(:),intent(in)::ramaseq
  !  integer,intent(in)::nres,nnode
  !  integer,dimension(nres),intent(in)::seq,resseq
  !  real,dimension(3,nres),intent(in)::angles,calpha
  !  real,dimension(nnode,nres),intent(in)::xgamma
  !  !type(model),intent(in)::modelR
  !  real,dimension(20,nres),intent(inout)::profile
  !  character::rama_c,aachar
  !  character(len=6)::sschar="HEGST_"
  !  character(len=mrama)::rama_string="HGBEdbeLlxc"
  !  character(len=1000)::ramafile
  !  character(len=300)::aline
  !  character(len=20),parameter::aa="ACDEFGHIKLMNPQRSTVWY"
  !  integer::i,ios,gunit=29,runit=30,ires,iaa,inode
  !  real,dimension(20)::gprofile
  !  real,dimension(mdssp)::ss
  !  real,dimension(0:mrama,20)::ramaprofile
  !  ramafile="ramaprofile.txt"
  !  ramaprofile=1./20.
  !  open(runit,file=trim(ramafile),status='old',form='formatted',iostat=ios)
  !  if(ios==0)then
  !    do
  !      read(runit,'(a)',iostat=ios)aline;if(ios/=0)exit;if(aline(1:1)=="!")cycle
  !      read(aline,*)rama_c
  !      do i=1,mrama
  !        if(rama_string(i:i)==rama_c)exit
  !      enddo
  !      if(i<mrama+1)read(aline,*)rama_c,ramaprofile(i,1:20)
  !    enddo
  !  endif
  !  close(runit)
  !  open(gunit,file=trim(gcafile),status='replace',form='formatted',iostat=ios)
  !  if(ios/=0)then
  !    write(*,*)"iostat=",ios
  !    stop 'hmmstrgca :: error opening output gca file'
  !  endif
  !  write(gunit,'(a)')'HMM: '//trim(mfile)//' prob_flag= 1 num_nodes= 282'
  !  !!GCA:  153l_ _ _ 185 Filename= /home/bystrc/hmm/db/val_combined.12mar99
  !  write(gunit,'(a,a,i5," Filename= ",a)')"GCA:  ",gcafile(1:4)//" _ _",nres,trim(gcafile)
  !  do ires=1,nres
  !    !!Write RESIDUE line
  !    aachar="X"
  !    if(seq(ires)>0.and.seq(ires)<=20)aachar=aa(seq(ires):seq(ires))
  !    write(gunit,'(a,2i4,1x,a1,$)')'RESIDUE: ',ires,resseq(ires),aachar
  !    write(gunit,'(3f9.1,3f9.2,a)')angles(1:3,ires),calpha(1:3,ires),' . . .'
  !    !!Write SS line
  !    ss(1:mdssp)=0.
  !    do inode=1,nnode
  !      ss(1:mdssp)=ss(1:mdssp)+xgamma(inode,ires)*modelR%d(inode,1:mdssp)
  !    enddo
  !    ss(1:mdssp)=ss(1:mdssp)/sum(ss(1:mdssp))
  !    iaa=maxloc(ss(1:mdssp),dim=1)
  !    aachar=sschar(iaa:iaa)
  !    write(gunit,'("SS: ",a1,6f9.6)')aachar,ss(1:mdssp)
  !    !!Write PROFILE line
  !    !profile=0.
  !    !profile(seq(ires))=1.0
  !    write(gunit,'("PROFILE:",20f9.6)')profile(1:20,ires)
  !    !!Write RAMA line
  !    write(gunit,'("RAMA: ",a1)')ramaseq(ires)
  !    !!Write RAMAPROFILE line
  !    write(gunit,'("RAMAPROFILE:",20f9.6)')ramaprofile(index(rama_string,ramaseq(ires)),1:20)
  !    !!Write GPROFILE line
  !    gprofile(1:20)=0.
  !    do inode=1,nnode
  !      gprofile(1:20)=gprofile(1:20)+xgamma(inode,ires)*modelR%b(inode,1:20)
  !    enddo
  !    gprofile(1:20)=gprofile(1:20)/sum(gprofile(1:20))
  !    write(gunit,'("GPROFILE:",20f9.6)')gprofile
  !    !!Write GAMMA line
  !    write(gunit,'("GAMMA:",999f9.6)')xgamma(1:nnode,ires)
  !  enddo
  !  close(gunit)
  !end subroutine write_gca
!-------------------------------------------------------------------------------------------------------------
  subroutine release_model()
   implicit none
   if(allocated(modelR%prior)) deallocate(modelR%prior)
   if(allocated(modelR%dih))   deallocate(modelR%dih)
   if(allocated(modelR%ss))    deallocate(modelR%ss)
   if(allocated(modelR%a))     deallocate(modelR%a)
   if(allocated(modelR%b))     deallocate(modelR%b)
   if(allocated(modelR%r))     deallocate(modelR%r)
   if(allocated(modelR%d))     deallocate(modelR%d)
   if(allocated(modelR%c))     deallocate(modelR%c)
   if(allocated(modelR%logb))  deallocate(modelR%logb)
  end subroutine release_model
!------------------------------------------------------------------------------------------------------------
  !subroutine set_flags(ramaseq,nres,use_str,angles,iflag)
  subroutine set_flags(ramaseq,nres,angles,iflag)
   !----------------------------------------
   implicit none
   !----------------------------------------
   character(len=1),dimension(:),intent(out)::ramaseq
   integer,intent(in)::nres
   integer,dimension(:),intent(out)::iflag
   !logical,intent(in)::use_str
   real,dimension(:,:),intent(in)::angles
   !----------------------------------------
   integer::i
   !----------------------------------------
   ramaseq="?"
   iflag=1
   !if(use_str) then
   if(.True.) then
     do i=1,nres
        if(all(angles(:,i)<998.)) then
           iflag(i)=3
        endif
        if(iflag(i)==3) then
           ramaseq(i)=ramatype(angles(1,i),angles(2,i),angles(3,i))
        endif
     enddo
   endif
  end subroutine set_flags
!------------------------------------------------------------------------------------------------------------
  !subroutine getbeta(ct,beta,seqres_len,modelR,profile,outtrans,ramaseq,iflag)
  subroutine getbeta(ct,beta,seqres_len,profile,ramaseq,iflag)
   !--------------------------------------------
   ! Calculate the backwards algorithm, return beta values
   ! for each sequence position, each Markov state.
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   integer,intent(in)::seqres_len
   integer,dimension(:),intent(in)::iflag
   real,dimension(:),intent(inout)::ct
   real,dimension(:,:),intent(out)::beta
   real,dimension(:,:),intent(in)::profile
   !type(model),intent(in)::modelR
   !type(trans),dimension(:),intent(in)::outtrans
   !--------------------------------------------
   integer::i
   integer::j
   integer::t
   integer::nq
   integer::jnd
   real(kind=8)::drat
   real(kind=8)::rrat
   real(kind=8)::ss
   real(kind=8)::cct
   !--------------------------------------------
   !!beta(nq,nres)  !!outtrans(nr)%c(nonzero)! profile(a(1:20),nres)
   ! Initialize all beta to zero
   !--------------------------------------------
   beta=0
   nq=modelR%n
   !--------------------------------------------
   ! Set all beta for the last position in the sequence
   ! equal to a constant. The constant is the scale factor ct
   ! which is output by the forward algorithm, getalpha.
   !--------------------------------------------
   do i=2,nq
     beta(i,seqres_len)=ct(seqres_len)
   enddo
   !--------------------------------------------
   ! State 1 is the "nought state". a non-emitting
   ! Markov state that connects other states. 
   ! Sum over all out-transitions from the nought state
   ! to every state jnd at the last position (seqres_len)
   ! to get beta for the naught state.
   !--------------------------------------------
   ss=0
   do j=1, outtrans(1)%c(0)
     jnd = outtrans(1)%c(j)
     ss = ss + modelR%a(1,jnd)*beta(jnd,seqres_len)* &
          !flagrat_brdc_prod(iflag(seqres_len),jnd,seqres_len,profile,modelR,ramaseq(seqres_len))
          flagrat_brdc_prod(iflag(seqres_len),jnd,seqres_len,profile,ramaseq(seqres_len))
   enddo
   beta(1,seqres_len)=ss
   !--------------------------------------------
   ! Stepping backward through the sequence positions
   !--------------------------------------------
   do t=seqres_len-1,1,-1
     cct=0
     !--------------------------------------------
     ! For each state i except the naught state
     ! sum the beta value. ct() are scale factors from
     ! the forward algorithm
     !--------------------------------------------
     do i=2,nq
       ss=0
       do j=1,outtrans(i)%c(0)
         jnd = outtrans(i)%c(j)
         ! bug fix, iflag(t) changed to iflag(t+1)
         !ss = ss + modelR%a(i,jnd)*beta(jnd,t+1)*flagrat_brdc_prod(iflag(t+1),jnd,t+1,profile,modelR,ramaseq(t+1))
         ss = ss + modelR%a(i,jnd)*beta(jnd,t+1)*flagrat_brdc_prod(iflag(t+1),jnd,t+1,profile,ramaseq(t+1))
       enddo
       beta(i,t)=ss*ct(t)
       cct = cct+beta(i,t)
     enddo
     !--------------------------------------------
     ! For the naught state, sum the beta value
     ! Note that the naught state transitions forward
     ! to t, not to t+1.
     !--------------------------------------------
     ss=0
     do j=1,outtrans(1)%c(0)
       jnd = outtrans(1)%c(j)
       !ss = ss + modelR%a(1,jnd)*beta(jnd,t)*flagrat_brdc_prod(iflag(t),jnd,t,profile,modelR,ramaseq(t))
       ss = ss + modelR%a(1,jnd)*beta(jnd,t)*flagrat_brdc_prod(iflag(t),jnd,t,profile,ramaseq(t))
     enddo
     beta(1,t)=ss
     if(cct==0)then
       write(*,*)"NO out going transistions!"
       stop
     endif
   enddo
  end subroutine
!-------------------------------------------------------------------------------------------------------------
  !subroutine getalpha(ct,alpha,seqres_len,modelR,profile,intrans,ramaseq,iflag)
  subroutine getalpha(ct,alpha,seqres_len,profile,ramaseq,iflag)
   !--------------------------------------------------
   ! This furnction calculates the forward algorithm,
   ! which is the probability of being in state q and
   ! position t in the sequence, given the HMM.
   !--------------------------------------------------
   implicit none
   !--------------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   integer,intent(in)::seqres_len
   integer,dimension(:),intent(in)::iflag
   real,dimension(:,:),intent(out)::alpha
   real,dimension(:,:),intent(in)::profile
   real,dimension(:),intent(inout)::ct
   !type(model),intent(in)::modelR
   !type(trans),dimension(:),intent(in)::intrans
   !--------------------------------------------------
   integer::i
   integer::j
   integer::t
   integer::nq
   integer::ind
   real(kind=8)::drat
   real(kind=8)::rrat
   real(kind=8)::ss
   real(kind=8)::junk
   real(kind=8)::tmp
   !--------------------------------------------------
   ! Initialize  forward values alpha as the
   ! probability of transition from the non-emitting
   ! "naught" state, which we use as the "begin" state.
   !--------------------------------------------------
   alpha=0;ct=0
   drat=1.;rrat=1.0
   ct(1)=0
   nq=modelR%n
   do j=1,nq
     !alpha(j,1)=modelR%a(1,j) * flagrat_brdc_prod(iflag(1),j,1,profile,modelR,ramaseq(1))
     alpha(j,1)=modelR%a(1,j) * flagrat_brdc_prod(iflag(1),j,1,profile,ramaseq(1))
     ct(1) = ct(1) + alpha(j,1)
   enddo
   !--------------------------------------------------
   ! SDcale factors ct are calculated such that the scaled
   ! alpha values are a probability distribution over
   ! nq states.
   !--------------------------------------------------
   ct(1) = 1./ct(1)
   do j=1,nq
     alpha(j,1) = alpha(j,1)*ct(1)
   enddo
   !--------------------------------------------------
   ! Calculate alpha for the naught state at t=1
   !--------------------------------------------------
   j=1
   ss=0
   do i=1, intrans(j)%c(0)
     ind = intrans(j)%c(i)
     ss = ss + alpha(ind,1)*modelR%a(ind,j)
   enddo
   alpha(j,1)=ss
   !--------------------------------------------------
   ! step through sequence positions
   !--------------------------------------------------
   do t=2,seqres_len
     ct(t)=0.
     !--------------------------------------------------
     ! calculate alpha for normal HMM states.
     !--------------------------------------------------
     do j=2,nq  
       ss=0
       do i=1, intrans(j)%c(0)
         ind = intrans(j)%c(i)
         ss = ss + alpha(ind,t-1)*modelR%a(ind,j)
       enddo
       !alpha(j,t) = ss * flagrat_brdc_prod(iflag(t),j,t,profile,modelR,ramaseq(t))
       alpha(j,t) = ss * flagrat_brdc_prod(iflag(t),j,t,profile,ramaseq(t))
       ct(t) = ct(t) + alpha(j,t)
     enddo
     if(ct(t)==0) then
       write(*,*) "No incoming transistions were found!"
       ct=1.
       return
     else
       ct(t) = 1./ct(t)
       do j=2,nq
         alpha(j,t) = alpha(j,t)*ct(t)
       enddo
       !--------------------------------------------------
       ! calculate alpha for naught state
       !--------------------------------------------------
       j=1
       ss=0
       do i=1,intrans(j)%c(0)
         ind = intrans(j)%c(i)
         ss = ss + alpha(ind,t) * modelR%a(ind,j)
       enddo
       alpha(j,t)=ss
     endif   
   enddo
  end subroutine
!-------------------------------------------------------------------------------------------------------------
  character(len=1) function ramatype(ph,ps,om)
   !--------------------------------------------
   ! This function assigns a character based on
   ! input phi/psi angles.
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   character(len=mrama)::rama_string="HGBEdbeLlxc"
   real,intent(in)::ph
   real,intent(in)::ps
   real,intent(in)::om
   !--------------------------------------------
   integer::i
   integer::icen
   real::d
   real::dmin
   real::ds
   real::df
   real,dimension(11,2)::ramacen
   !--------------------------------------------
   ! ramacen(:,1:2) are phi, psi angles for
   ! centroids of vononoi spaces in the
   ! Ramachandran plot. Voronoi spaces
   ! are calculated using Pythagoean distance
   ! in toroidal angle space.
   !--------------------------------------------
   ramacen(1,1)=-61.91   ;ramacen(1,2)=-45.20
   ramacen(2,1)=-109.78  ;ramacen(2,2)=20.88
   ramacen(3,1)=-70.58   ;ramacen(3,2)=147.22
   ramacen(4,1)=-132.89  ;ramacen(4,2)=142.43
   ramacen(5,1)=-135.03  ;ramacen(5,2)=77.26
   ramacen(6,1)=-85.03   ;ramacen(6,2)=72.26
   ramacen(7,1)=-165.00  ;ramacen(7,2)=175.00
   ramacen(8,1)=55.88    ;ramacen(8,2)=38.62
   ramacen(9,1)=85.82    ;ramacen(9,2)=-.03
   ramacen(10,1)=80.     ;ramacen(10,2)=-170.00
   ramacen(11,1)=-70.0   ;ramacen(11,2)=150.00
   if(ph==999.)then; ramatype="?"; return; endif
   if(ps==999.)then; ramatype="?"; return; endif
   if(om<90.and.om>-90)then; ramatype="c"; return; endif
   dmin=999
   do i=1,mrama-1
     df=abs(ph-ramacen(i,1)); if(df>180.) df= 360.-df  
     ds=abs(ps-ramacen(i,2)); if(ds>180.) ds= 360.-ds  
     d=sqrt(ds*ds+df*df)
     if(d<dmin) then
       icen=i; dmin=d
     endif  
   enddo
   ramatype=rama_string(icen:icen)
  end function ramatype
!-------------------------------------------------------------------------------------------------------------
  !real function bprof_bg(iq,t,profile,modelR)
  real function bprof_bg(iq,t,profile)
   !---------------------------
   ! This function calculates the profile*profile
   ! score as p log(q/b) where p is the sequence
   ! profile and q is the state profile and b is
   ! the background profile.
   ! A constant amount EPS of background profile
   ! is mixed in to prevent log(0).
   ! The value returned is a probability ratio.
   !---------------------------
   implicit none
   !---------------------------
   integer,intent(in)::t
   integer,intent(in)::iq
   real,dimension(:,:),intent(in)::profile
   !type(model),intent(in)::modelR
   !---------------------------
   integer::a
   real(kind=8)::logsum
   real(kind=8)::tmp
   !---------------------------
   logsum=0
   do a=1,20
     if(profile(a,t)/=0)then
       tmp = log((modelR%b(iq,a)+(EPS*modelR%bground(a)))/((1.+EPS)*modelR%bground(a)))    
       logsum = logsum + profile(a,t)*ncount*tmp   
     endif
   enddo
   bprof_bg=exp(logsum)
  end function bprof_bg
!-------------------------------------------------------------------------------------------------------------
  !real function flagrat_brdc_prod(jflag,iq,t,profile,modelR,rama_c)
  real function flagrat_brdc_prod(jflag,iq,t,profile,rama_c)
   !-----------------------------------------
   ! This function returns the unconditional proability of state iq
   ! at position t given the sequence+structure data (profile,
   ! rama_c) jflag determines which data is used.
   ! Returned value is a joint probability ratio.
   !-----------------------------------------
   implicit none
   !-----------------------------------------
   character(len=1),intent(in)::rama_c
   integer,intent(in)::jflag
   integer,intent(in)::t,iq
   real,dimension(:,:),intent(in)::profile
   !type(model),intent(in)::modelR
   !-----------------------------------------
   integer::k
   real::brat
   real::rrat
   character(len=mrama)::rama_string="HGBEdbeLlxc"
   !-----------------------------------------
   brat=1; rrat=1
   if(iq==1) then
     flagrat_brdc_prod=brat
     return
   endif
   if(jflag==3)then
     do k=1,mrama
       if(rama_string(k:k)==rama_c) exit
     enddo
     if( k < mrama+1) then
       rrat = (modelR%r(iq,k) + EPS*modelR%rground(k))/((1.+EPS)*modelR%rground(k))
     !else
     !  write(*,*)"RAMA character not fount! ",rama_c
     !  stop
     endif
   endif  
   if(jflag<=3)then
     brat = bprof_bg(iq,t,profile)
   endif
   flagrat_brdc_prod = brat*rrat
  end function flagrat_brdc_prod
!-------------------------------------------------------------------------------------------------------------
  !subroutine get_intrans(intrans,modelR)
  subroutine get_intrans()
   !---------------------------------------
   implicit none
   !---------------------------------------
   !type(model),intent(in)::modelR
   !type(trans),dimension(:),allocatable,intent(out)::intrans
   !---------------------------------------
   integer::nq
   integer::ncount
   integer::i
   integer::j
   integer::ios
   !---------------------------------------
   nq=modelR%N
   allocate(intrans(nq))
   do j=1,nq
     ncount=0 
     do i=1,nq
       if(modelR%a(i,j)/=0) ncount = ncount + 1
     enddo
     allocate(intrans(j)%c(0:ncount))
     ncount=0 
     do i=1,nq
       if(modelR%a(i,j)/=0) then
         ncount = ncount + 1
         intrans(j)%c(ncount) = i
       endif
       intrans(j)%c(0) = ncount
     enddo  
   enddo
  end subroutine get_intrans
!-------------------------------------------------------------------------------------------------------------
  !subroutine get_outtrans(outtrans,modelR)
  subroutine get_outtrans()
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   !type(model),intent(in)::modelR
   !type(trans),dimension(:),allocatable,intent(out)::outtrans
   !--------------------------------------------
   integer::nq
   integer::ncount
   integer::ios
   integer::i
   integer::j
   !--------------------------------------------
   nq=modelR%N
   allocate(outtrans(nq))
   do i=1,nq
     ncount=0 
     do j=1,nq
       if(modelR%a(i,j)/=0) ncount = ncount + 1
     enddo
     allocate(outtrans(i)%c(0:ncount))
     ncount=0 
     do j=1,nq
       if(modelR%a(i,j)/=0) then
         ncount = ncount + 1
         outtrans(i)%c(ncount) = j
       endif
       outtrans(i)%c(0) = ncount
     enddo  
   enddo
  end subroutine get_outtrans
!-------------------------------------------------------------------------------------------------------------
  !subroutine readmodel(mfile,modelR,sumd,sumr,sumc)
  !subroutine readmodel(mfile,sumd,sumr,sumc)
  subroutine readmodel(mfile,sumd,sumr,sumc)
   !---------------------------------------
   implicit none
   !---------------------------------------
   character(len=1000000),intent(in)::mfile
   real,dimension(:),intent(out)::sumd
   real,dimension(:),intent(out)::sumr
   real,dimension(:),intent(out)::sumc
   !type(model),intent(out)::modelR
   !---------------------------------------
   character(len=300)::aline
   character(len=200)::junk
   integer::ios
   integer::i
   integer::j
   integer::nnodes
   integer::n
   real::x
   real::z
   !---------------------------------------
 
   open(1,file=mfile,status='old',iostat=ios)
   if(ios/=0) then 
      write(*,*) "Cannot open model file!"
      stop
   endif
   sumr=0; sumd=0; sumc=0
   do 
     read(1,'(a)',iostat=ios) aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm","1"
       stop
     endif 
     if(aline(1:16) == "default_bg_freqs")then
        read(aline(17:197),*)modelR%bground(1:20)
        write(*,'(a,20f9.5)')"BACKGROUND FREQS ",modelR%bground(1:20)
     endif
     if(aline(1:9) == "num_nodes")then
        read(aline(10:15),*)modelR%N
        write(*,'(a,i5)')"NUM NODES ",modelR%N
        exit
     endif
   enddo
   modelR%N = modelR%N + 1
   n=modelR%N
   allocate(modelR%prior(n))
   allocate(modelR%dih(n,3))
   allocate(modelR%ss(n))
   allocate(modelR%a(n,n))
   allocate(modelR%b(n,m))
   allocate(modelR%r(n,mrama))
   allocate(modelR%d(n,mdssp))
   allocate(modelR%c(n,mctxt))
   allocate(modelR%logb(n,m))
   do 
     read(1,'(a)',iostat=ios) aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm","2"
       stop
     endif 
     if(aline(1:9) == "unk_node ")then
        read(aline(14:207),*)modelR%prior(1)
        read(1,*)modelR%b(1,1:m)
        modelR%dih(1,1)= -75
        modelR%dih(1,2)= -15
        modelR%dih(1,3)= 180
        modelR%ss(1) = "_"
     endif
     if(aline(1:11) == "unk_node_ss") then
        read(aline(19:63),*)modelR%d(1,1:mdssp)
        sumd = sumd + modelR%prior(1)*modelR%d(1,:)
     endif
     if(aline(1:13) == "unk_node_rama") then
        read(aline(21:105),*)modelR%r(1,1:mrama)
        sumr = sumr + modelR%prior(1)*modelR%r(1,:)
        exit
     endif            
   enddo
   !!! For all known nodes
   do i=2,n
     read(1,'(a)',iostat=ios)aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm","3"
       stop
     endif 
     read(aline(6:8),*)z
     if(i-1/=z .or. aline(1:4) /= "node") then
       write(*,*) "Improperly formatted modelR.hmm","4"
       stop
     endif 
     read(aline,*)junk,junk,junk,junk,modelR%prior(i),modelR%dih(i,1:3)
     read(1,*)modelR%b(i,1:10)
     read(1,*)modelR%b(i,11:20)
     read(1,'(a)',iostat=ios)aline
     read(aline(21:66),*)modelR%d(i,1:mdssp)
     sumd = sumd + modelR%prior(i)*modelR%d(i,:)
     read(1,'(a)',iostat=ios)aline
     read(aline(24:108),*)modelR%r(i,1:mrama)
     sumr = sumr + modelR%prior(i)*modelR%r(i,:)
     read(1,'(a)',iostat=ios)aline
     read(aline(27:103),*)modelR%c(i,1:mctxt)
     sumc = sumc + modelR%prior(i)*modelR%c(i,:)
   enddo
   read(1,*)aline
   if(aline(1:13)/="transit_freqs") then
     write(*,*) "Improperly Formatted modelR.hmm"
   endif
   x=0; x = sum(sumr)
   modelR%rground = sumr/x
   x=0; x = sum(sumd)
   modelR%dground = sumd/x
   x=0; x = sum(sumc)
   modelR%cground = sumc/x
   do i=1,modelR%N
     read(1,*)modelR%a(i,1:modelR%n)
     !write(*,*)modelR%a(i,1:n)
   enddo
   close(1)
  end subroutine readmodel
!-------------------------------------------------------------------------------------------------------------
!  subroutine read_seq(sfile,seq,seqnres)
!   !-------------------------------------
!   implicit none
!   !-------------------------------------
!   integer,intent(out)::seqnres
!   character(len=*),intent(in)::sfile
!   character(len=*),intent(out)::seq
!   !-------------------------------------
!   character(len=50)::aline
!   integer::ios
!   !-------------------------------------
!   open(33,file=sfile,iostat=ios)
!   if(ios/=0)then
!     write(*,*)"Cannot open seq file!"
!     stop
!   endif
!   read(33,*)aline
!   if(aline(1:1)/=">")then
!      write(*,*)"Sequence file should be fasta format"
!      stop
!   endif
!   read(33,*)aline
!   seq = trim(aline)
!   do 
!     read(33,*,iostat=ios)aline
!     if(ios/=0)exit
!     seq = trim(seq)//trim(aline)
!   enddo
!   seqnres=len_trim(seq)
!   close(33)
!  end subroutine read_seq
!-------------------------------------------------------------------------------------------------------------
  subroutine read_seq(sfile,seq,seqnres)
    implicit none
    integer,intent(out)::seqnres
    character(len=*),intent(in)::sfile
    character(len=3000)::seqch
    integer,dimension(3000),intent(out) ::seq
    character(len=200)::aline
    character(len=20),parameter::aa1="ACDEFGHIKLMNPQRSTVWY"
    integer::ios,ires
    open(33,file=sfile,iostat=ios)
    if(ios/=0)then
      write(*,*)"Cannot open seq file"
      stop
    endif
    read(33,*)aline
    if(aline(1:1)/=">")then
      write(*,*)"Sequence file should be fasta format"
      stop
    endif
    read(33,*)aline
    seqch=trim(aline)
    do
      read(33,*,iostat=ios)aline
      if(ios/=0)exit
      seqch=trim(seqch)//trim(aline)
    enddo
    seqnres=len_trim(seqch)
    do ires=1,seqnres
      seq(ires)=index(aa1,seqch(ires:ires))
    enddo
    close(33)
  end subroutine read_seq
!-------------------------------------------------------------------------------------------------------------
  subroutine read_profile(pfile, profile, seq, nres)
    implicit none
    integer, intent(in)::nres
    character(len=*), intent(in)::pfile
    integer, dimension(3000), intent(in)::seq
    real, dimension(20, nres), intent(inout)::profile
    integer::i, j, ios, ires
    character(len=300)::aline
    character(len=3000)::seqch
    real, dimension(nres)::tmpvec
    character(len=1), dimension(21)::res1=(/'A','C','D','E','F','G','H','I','K','L',&
                                            'M','N','P','Q','R','S','T','V','W','Y',&
                                            'X'/)
    character(len=20), parameter::aa1="ACDEFGHIKLMNPQRSTVWY"

    do ires=1, nres
      !! Change any trailing positions (waters) to 21's.
      !! CB Tue Dec 28 15:17:31 EST 2010
      i=seq(ires)
      if (i<=0) i=21
      seqch(ires:ires)=res1(i)
    end do
    open (1, file=pfile, form="formatted", status="old", iostat=ios)
    if (ios/=0) then
      write (*, *) "subroutine read_profile (hmmstr.f90) :: cannot open profile"
      stop
    end if
    read (1, '(A)') aline
    do while (aline(12:21)/="A  R  N  D")
      read (1,'(A)',iostat=ios) aline
	if(ios/=0)exit
    end do
    profile=0.0
    do i=1, nres
      read (1, '(A)', iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:5)=="     ") exit
      if (aline(7:7)/=seqch(i:i)) then
        write (*, *) i, aline(7:7), seqch(i:i)
        write (*, *) "subroutine read_profile (hmmstr.f90) :: profile and seqfile nres do not match"
        stop
      end if
      read (aline(8:), *) profile(1:20, i), profile(1:20, i)
    end do
    if (i>nres+1) then
      write (*, *) "subroutine read_profile (hmmstr.f90) :: profile and seqfile nres do not match"
      stop
    end if
    close(1)
    !!THIS IS A CRAZY WAY to permute the AA order.
    tmpvec=profile(2, :); profile(2, :)=profile(5, :); profile(5, :)=tmpvec
    tmpvec=profile(3, :); profile(3, :)=profile(4, :); profile(4, :)=tmpvec
    tmpvec=profile(4, :); profile(4, :)=profile(7, :); profile(7, :)=tmpvec
    tmpvec=profile(5, :); profile(5, :)=profile(14, :); profile(14, :)=tmpvec
    tmpvec=profile(6, :); profile(6, :)=profile(8, :); profile(8, :)=tmpvec
    tmpvec=profile(7, :); profile(7, :)=profile(9, :); profile(9, :)=tmpvec
    tmpvec=profile(8, :); profile(8, :)=profile(10, :); profile(10, :)=tmpvec
    tmpvec=profile(9, :); profile(9, :)=profile(12, :); profile(12, :)=tmpvec
    tmpvec=profile(10, :); profile(10, :)=profile(11, :); profile(11, :)=tmpvec
    tmpvec=profile(11, :); profile(11, :)=profile(13, :); profile(13, :)=tmpvec
    tmpvec=profile(13, :); profile(13, :)=profile(15, :); profile(15, :)=tmpvec
    tmpvec=profile(14, :); profile(14, :)=profile(15, :); profile(15, :)=tmpvec
    tmpvec=profile(18, :); profile(18, :)=profile(19, :); profile(19, :)=tmpvec
    tmpvec=profile(20, :); profile(20, :)=profile(18, :); profile(18, :)=tmpvec
    do i=1, nres
      if (all(profile(1:20, i)==0)) then
        do j=1, 21
          if (res1(j)==seqch(i:i)) exit
        end do
        if (j/=21) then
          profile(j, i)=1.0
        else
          profile(1:20, i)=modelR%bground(1:20)
        end if
      end if
      profile(1:20, i)=profile(1:20, i)/sum(profile(1:20, i))
      !write (*, *) profile(1:20, i)
    end do
  end subroutine read_profile
!-------------------------------------------------------------------------------------------------------------
  subroutine get_angles(backbone,nres,angles)
   !---------------------------------------------------
   implicit none
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::backbone
   real,dimension(:,:),intent(out)::angles
   !---------------------------------------------------
   integer::i
   integer::j
   integer::iatom
   real,dimension(3,4)::phixyz
   real,dimension(3,4)::psixyz
   real,dimension(3,4)::omgxyz
   !---------------------------------------------------
   angles=999.
   i=1
   psixyz(:,1:3) = backbone(1:3,i:i+2) 
   psixyz(:,4) = backbone(1:3,i+4) 
   omgxyz(:,1:2) = backbone(1:3,i+1:i+2) 
   omgxyz(:,3:4) = backbone(1:3,i+4:i+5) 
   angles(2,1) = cal_tors(psixyz)
   angles(3,1) = cal_tors(omgxyz)
   j=2
   do i=5,(nres-1)*4,4
     phixyz(:,1) = backbone(:,i-2) 
     phixyz(:,2:4) = backbone(:,i:i+2) 
     psixyz(:,1:3) = backbone(1:3,i:i+2) 
     psixyz(:,4) = backbone(1:3,i+4) 
     omgxyz(:,1:2) = backbone(1:3,i+1:i+2) 
     omgxyz(:,3:4) = backbone(1:3,i+4:i+5) 
     angles(1,j) = cal_tors(phixyz)
     angles(2,j) = cal_tors(psixyz)
     angles(3,j) = cal_tors(omgxyz)
     j = j + 1
   enddo
   phixyz(:,1) = backbone(:,i-2) 
   phixyz(:,2:4) = backbone(:,i:i+2)
   angles(1,j) = cal_tors(phixyz)
  end subroutine get_angles
!-------------------------------------------------------------------------------------------------------------
  real function cal_tors(xyz) 
   !----------------------------------
   implicit none
   !----------------------------------
   real,dimension(3,4),intent(in)::xyz
   !----------------------------------
   integer::i
   real::angle
   real::handcheck
   real::len_c13
   real::len_c24
   real,dimension(3)::vec
   real,dimension(3)::c13
   real,dimension(3)::c24
   real,dimension(3)::c13xc24
   real,dimension(3)::rotaxis
   real,dimension(3,4)::tmpxyz
   real,parameter::pi=3.14159265
   real,parameter::rad2deg=180/pi
   !----------------------------------
   vec=xyz(:,2)
   do i=1,4
     tmpxyz(:,i) = xyz(:,i) - vec
   enddo
   call crossprod(tmpxyz(:,1),tmpxyz(:,3),c13)
   rotaxis=tmpxyz(:,3)
   vec=xyz(:,3)
   do i=1,4
     tmpxyz(:,i) = xyz(:,i) - vec
   enddo
   call crossprod(tmpxyz(:,2),tmpxyz(:,4),c24)
   call crossprod(c24,c13,c13xc24)
   handcheck = sum(c13xc24*rotaxis)
   len_c13 = sqrt(sum(c13*c13));  len_c24 = sqrt(sum(c24*c24))
   angle = sum(c13*c24)/(len_c13*len_c24)
   if(angle>1) angle=1;   if(angle < -1) angle = -1 !! Safety first.
   angle = acos(angle)
   if(handcheck>0) then
     angle = -angle 
   endif
   cal_tors = angle*rad2deg
  end function cal_tors
!-------------------------------------------------------------------------------------------------------------
  !subroutine read_backbone(pdbfile,chain,bbxyz,seq,nres,resseq)
  subroutine read_backbone(pdbfile,chain,seq,nres,resseq)
    implicit none
    character(len=*),intent(in)::pdbfile
    integer,dimension(nres),intent(out)::seq
    integer,dimension(nres),intent(out)::resseq
    character(len=3000)::seqch
    character(len=1),intent(inout)::chain
    integer,intent(inout)::nres
    !type(bb),dimension(nres),intent(inout)::bbxyz
    character(len=1)::altloc
    character(len=100)::aline
    integer::i
    integer::j
    integer::k
    integer::iatom
    integer::ios
    integer::ires,mres
    integer::atype
    character(len=6)::last
    character(len=9)::curres
    character(len=9)::qres
    character(len=3),dimension(20)::three=&
      (/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',&
        'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'/)
    character(len=1),dimension(21)::res1=&
      (/'A','C','D','E','F','G','H','I','K','L',&
        'M','N','P','Q','R','S','T','V','W','Y','X'/)
    open(2,file=pdbfile,status='old',iostat=ios)
    if(ios/=0)then
      write(*,*)"Cannot open pdbfile ",trim(pdbfile)
      stop
    endif
    ires=0;curres=" ";seqch=" ";seq=0
    do
      read(2,'(a)',iostat=ios)aline
      if(ios/=0)exit
      if(aline(1:6)=='ENDMDL')exit                            !!Use first model of NMR structures
      if(aline(1:3)=='TER'.and.ires/=0)exit                   !!Atoms with the correct chain id but not part of the chain
      if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM')cycle
      if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
      if(aline(18:20) == 'HOH')cycle
      if(aline(22:22)==" ")chain="_"
      !if(aline(13:16)/=' CA ')cycle                          !!Only one CA per residue allowed
      if(aline(18:26)==curres)cycle
      curres=aline(18:26)
      ires=ires+1
      read(aline(23:26),*)resseq(ires)
      k=21
      do j=1,20
        if(aline(18:20)==three(j))then
          k=j;exit
        endif
      enddo
      seqch(ires:ires)=res1(k)
      seq(ires)=k
    enddo
    mres=ires;iatom=0;ires=0
    do ires=1,nres
      do iatom=1,4
        bbxyz(ires)%xyz(1:3,iatom)=999.
      enddo
    enddo
    rewind(2)
    curres=" ";altloc=" ";ires=0;
    do
      iatom=0
      read(2,'(a)',iostat=ios)aline
      if(ios/=0)exit
      if(aline(1:6)=='ENDMDL')exit
      if(aline(1:3)=='TER'.and.ires/=0)exit                   !!Atoms with the correct chain id but not part of the chain
      if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM')cycle
      if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
      altloc=aline(17:17)
      if(aline(18:20) == 'HOH')cycle
      if(altloc/=" ".and.altloc/="A".and.altloc/="1")cycle
      if(aline(18:26)/=curres)then
        ires=ires+1                                           !!Attempting to read this residue.
        curres=aline(18:26)                                   !!Currently reading...
      endif
      atype=idatom(aline(13:16),aline(20:22))
      if(atype==-1)cycle                                      !!This is not a backbone atom. Reject.
      if(ires>nres)then
        write(*,*)"MISSING BACKBONE ATOMS?(1) ",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
        stop
      endif
      read(aline(31:54),'(3f8.3))')bbxyz(ires)%xyz(1:3,atype)
      iatom=iatom+1
      do
        if(iatom==4)exit                                      !!All four backbone atoms have been read
        read(2,'(a)',iostat=ios)aline
        if(ios/=0)exit
        if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM')cycle
        if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
        altloc=aline(17:17)
        if(altloc/=" ".and.altloc/="A".and.altloc/="1")cycle
        qres=aline(18:26)
        if(aline(18:20) == 'HOH')cycle                          !!Couldn't find all four atoms.
        if(curres/=qres)then
          write(*,*)"MISSING BACKBONE ATOMS?(2) ",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
          iatom=0
          ires=ires+1
          curres=aline(18:26)
        endif
        atype=idatom(aline(13:16),aline(20:22))
        if(atype==-1)cycle
        if(all(bbxyz(ires)%xyz(1:3,atype)>998.))then
          read(aline(31:54),'(3f8.3))')bbxyz(ires)%xyz(1:3,atype)
          iatom=iatom+1
        else
          write(*,*)"ATTEMPTING TO OVERWRITE ATOM COORDINATES ",ires," ",aline(18:26)," ",trim(pdbfile)," ",chain
          stop
        endif
      enddo
    enddo
    close(2)
  end subroutine read_backbone
!-------------------------------------------------------------------------------------------------------------
  subroutine crossprod(v1,v2,v3)
   !------------------------------
   implicit none
   !------------------------------
   real,intent(in)::v1(3)
   real,intent(in)::v2(3)
   real,intent(out)::v3(3)
   !------------------------------
   v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
   v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
   v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end subroutine crossprod
!-------------------------------------------------------------------------------------------------------------
  integer function idatom(atom,res)
   !----------------------------------
   implicit none
   !----------------------------------
   character(len=4),intent(in)::atom
   character(len=3),intent(in)::res
   !----------------------------------
   integer::i
   character(len=4),dimension(1:4)::atype=(/" N  "," CA "," C  "," O  "/)
   character(len=3),dimension(20)::three=(/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE',&
     'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL',&
     'TRP','TYR'/)
   character(len=1),dimension(21)::res1=(/'A','C','D','E','F','G','H','I','K','L','M','N',&
     'P','Q','R','S','T','V','W','Y','X'/)
   !----------------------------------
   do i=1,4
     if(atom==atype(i))exit
   enddo
   if(i==5)then
     idatom = -1
   else
     idatom=i
   endif
  end function idatom
!-------------------------------------------------------------------------------------------------------------
end module hmmstr

