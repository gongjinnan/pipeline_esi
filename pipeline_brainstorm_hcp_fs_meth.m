MriFile = 'H:\chang_xin_20170831_5453_010_MITCHELL_JORGE_3D_T1.nii_out';
addpath E:\fieldtrip-20170224
load sa_template;
Brainstorm_route = cd;

nVertices = 10000;
load('4_105923_MEG_3-Restin_rmegpreproc_elec.mat')
nVertices = 10000;

my_addjava(Brainstorm_route);

[sMri,sHead, Cortex ] = my_import_anatomy_fs( MriFile, nVertices, 0, [], 0);
vol.bnd.pos = Cortex.Vertices;
vol.bnd.tri = Cortex.Faces;
 mri = ft_read_mri('H:\chang_xin_20170831_5453_010_MITCHELL_JORGE_3D_T1.nii_out\mri\T1.nii');
 mri.coordsys = 'nifti';
sa_meg1=my_mk_sa_meg_mri(sa_template,vol,mri);

% checkmritrafo(sa_meg1,mri);

ori = elec_grad.chanori(strcmp('meg',elec_grad.chantype),:);

ref = [elec_grad.chanpos(strcmp('refmag',elec_grad.chantype),:), elec_grad.chanori(strcmp('refmag',elec_grad.chantype),:) ];


sens1 = [elec_pos,ori];
head.vc = sHead.Vertices;
head.tri = sHead.Faces;
figure;showsurface(head)
hold on;
show_megsystem(sens1,0.01);

xmin=min(elec_pos)-0.01;
xmax=max(elec_pos)+0.01;
axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)]); 

sa=my_mk_sa_meg_forward(sa_meg1,sens1,refs);
%  checkmritrafo(sa,mri)
%  checksenslocs(sa)
 L=grid2L(sa.grid_medium_indi,sa.fp_indi);
% LeadField=grid2L(Cortex.Vertices,fp);
% [ne,nv,dim] = size(LeadField);
% L = reshape(permute(LeadField,[1,3,2]), [ne,nv*dim]);
% L =  bst_gain_orient(L, Cortex.VertNormals);

