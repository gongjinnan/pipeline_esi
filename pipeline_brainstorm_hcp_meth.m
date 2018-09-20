MriFile = 'E:\workshop_beijing_2\anatomy\T1w_acpc_dc_restore.nii.gz';
FileFormatMri = 'ALL-MNI';
TransFile = 'E:\workshop_beijing_2\anatomy\105923_MEG_anatomy_transform.txt';
TessFileL = 'E:\workshop_beijing_2\anatomy\105923.L.midthickness.4k_fs_LR.surf.gii';
FileFormatTess = 'GII-MNI';
Brainstorm_route = 'E:\Brainstorm 2.0\brainstorm3';
TessFileR = 'E:\workshop_beijing_2\anatomy\105923.R.midthickness.4k_fs_LR.surf.gii';
DataFile = 'E:\workshop_beijing_2\100307\unprocessed\MEG\3-Restin\4D\c,rfDC';
nVertices = 10000;
erodeFactor = 0;
fillFactor =2;
[ sMri, Cortex, sHead ] = load_pipeline_brainstorm(MriFile,FileFormatMri,TransFile,TessFileL,TessFileR,FileFormatTess, nVertices,erodeFactor,fillFactor);


%BEM
nvert =[1922 1922 1922];
thickness= [7 4 3];
[ sInner,sOuter,sHeadBEM, sHead ] = bem_surfaces_brainstorm( sMri,sHead,Cortex,nvert,thickness,Brainstorm_route );

ChannelMat = my_in_fopen_4d(DataFile);
iMeg = good_channel(ChannelMat.Channel,[],'MEG');
elect = [ChannelMat.Channel(iMeg).Loc];
elect =elect(:,1:4:end)';

ori = [ChannelMat.Channel(iMeg).Orient];
ori =ori(:,1:4:end)';
sens1 = [elect,ori];
head.vc = sHead.Vertices;
head.tri = sHead.Faces;
figure;showsurface(head)
hold on;
show_megsystem(sens1,0.01);

xmin=min(elect)-0.01;
xmax=max(elect)+0.01;
axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)]); 

p1=15;
p2=10;

[vc,center,radius,coeffs]=pointsonsurface(Cortex.Vertices(:,1:3),Cortex.Vertices(:,1:3),p1);
fp=meg_ini(vc,center',p2,sens1);  
LeadField=grid2L(Cortex.Vertices,fp);
[ne,nv,dim] = size(LeadField);
L = reshape(permute(LeadField,[1,3,2]), [ne,nv*dim]);
L =  bst_gain_orient(L, Cortex.VertNormals);

