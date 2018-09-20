MriFile = 'H:\Pedro Ariel\MC0000003_t13d_anatVOL_20060108014323_2.nii_out';
Brainstorm_route = cd;


load('elec_pos_ICBM.mat')
my_addjava(Brainstorm_route);

nVertices = 10000;

[sMri, Cortex,WhiteLow, MidLow, sHead] = my_import_anatomy_fs( MriFile, nVertices, 0, [], 0);

%BEM
nvert =[1922 1922 1922];
thickness= [7 4 3];
[ sInner,sOuter,sHeadBEM, sHead ] = bem_surfaces_brainstorm(sMri,sHead,Cortex,nvert,thickness,Brainstorm_route);
[locsx] = channel_project_scalp(sHeadBEM.Vertices,elect);

iVertInside = find(inpolyhd(locsx', sInner.Vertices, sInner.Faces));
sa_in.vc = cell(1,3);
sa_in.vc{1}.vc = sInner.Vertices;
sa_in.vc{2}.vc = sOuter.Vertices;
sa_in.vc{3}.vc = sHeadBEM.Vertices;
sa_in.Cortex = Cortex.Vertices;
sa_in.head.vc = sHeadBEM.Vertices;
sa_in.head.tri = sHeadBEM.Faces;
LeadField=mk_sa_eeg_new(sa_in,locsx);

[ne,nv,dim] = size(LeadField);
L = reshape(permute(LeadField,[1,3,2]), [ne,nv*dim]);
L =  bst_gain_orient(L, Cortex.VertNormals);


