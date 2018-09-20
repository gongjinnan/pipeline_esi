MriFile = 'E:\FS_folder';
DataFile = 'E:\workshop_beijing_2\100307\unprocessed\MEG\3-Restin\4D\c,rfDC';
% TransFile = 'E:\workshop_beijing_2\anatomy\105923_MEG_anatomy_transform.txt';
% FileFormatMri = 'ALL-MNI';
% TessFileL = 'E:\workshop_beijing_2\anatomy\105923.L.midthickness.4k_fs_LR.surf.gii';
% FileFormatTess = 'GII-MNI';
Brainstorm_route = cd;
my_addjava(Brainstorm_route);
% TessFileR = 'E:\workshop_beijing_2\anatomy\105923.R.midthickness.4k_fs_LR.surf.gii';

nVertices = 10000;
erodeFactor = 0;
fillFactor =2;
%[ sMri, Cortex, sHead ] = load_pipeline_brainstorm(MriFile,FileFormatMri,TransFile,TessFileL,TessFileR,FileFormatTess, nVertices,erodeFactor,fillFactor);
[sMri, Cortex,WhiteLow, MidLow, sHead] = my_import_anatomy_fs( MriFile, nVertices, 0, [], 0);

%BEM
nvert =[1922 1922 1922];
thickness= [7 4 3];
[ sInner,sOuter,sHeadBEM, sHead ] = bem_surfaces_brainstorm( sMri,sHead,Cortex,nvert,thickness,Brainstorm_route );

ChannelMat = my_in_fopen_4d(DataFile);
elect = [ChannelMat.Channel(:).Loc];
elect =elect(:,1:4:end)';

ori = [ChannelMat.Channel(:).Orient];
ori =ori(:,1:4:end)';

head.vc = sHead.Vertices;
head.tri = sHead.Faces;
figure;showsurface(head)
hold on;
show_megsystem([elect,ori],0.01);

xmin=min(elect)-0.01;
xmax=max(elect)+0.01;
axis([xmin(1) xmax(1) xmin(2) xmax(2) xmin(3) xmax(3)]); 


OPTIONS.GridLoc = Cortex.Vertices;
OPTIONS.isEeg = 0;
OPTIONS.isMeg = 1;
OPTIONS.BemSurf = cell(3,1);
OPTIONS.BemSurf{3} = sInner;
OPTIONS.BemSurf{2}= sOuter;
OPTIONS.BemSurf{1} = sHeadBEM;

OPTIONS.BemNames = {'Scalp'    'Skull'    'Brain'};
OPTIONS.BemCond = [1.0000    0.0125    1.0000];

OPTIONS.isAdjoint = 1;
OPTIONS.isAdaptative = not(OPTIONS.isAdjoint);
OPTIONS.GridOrient = Cortex.VertNormals;

OPTIONS.Channel = ChannelMat.Channel;
OPTIONS.iMeg = good_channel(ChannelMat.Channel,[],'MEG');
iBad = find(any(isnan(OPTIONS.GridOrient),2) | any(isinf(OPTIONS.GridOrient),2) | (sqrt(sum(OPTIONS.GridOrient.^2,2)) < eps));
if ~isempty(iBad)
    OPTIONS.GridOrient(iBad,:) = repmat([1 0 0], length(iBad), 1);
end
Gain = my_bst_openmeeg(OPTIONS);

Gain =  bst_gain_orient(Gain, OPTIONS.GridOrient);
