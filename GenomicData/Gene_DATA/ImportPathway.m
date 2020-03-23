%clear all; close all
learn_downstream = 1;
learn_upstream = 0;

ipmode = input('Enter 1:Import Pathway   2:Generate Adjance Matrix    3:Show Network Statistics    4:Plot Node:   5: Show Parents/Children:  ');

import_pathway = 0;    generate_adjc_matrix = 2;     show_data = 3;    browse_node = 4;  show_network=5;
show_all_parents = 0;  %for browse_node
initvars
%%
if ipmode == import_pathway 
    load (GeneNameListFile, 'GeneDataBase');
    TXT = textread([cd, '\Pathway\pid_120912_pathway.tab.tsv'], '%s', 'delimiter', '\n', 'whitespace','');
    TXT = strtrim(TXT);
    nItems  = length(TXT);

    if ~isempty(dir([cd, '\Pathway\pathway.mat']))
        load([cd, '\Pathway\pathway.mat'], 'TXT',  'node', 'edge', 'nEdges', 'nNodes');
        if nItems - nEdges - nNodes <5
            fprintf('File is already imported, Are you sure to overwrite? Press a key ....'); pause; pause; pause; pause; 
        end
    end


    %*************************************%
    %   Find Nodes                        %
    %*************************************%
    j=0;
    for i = 1: nItems
        pathitem = strsplit(TXT{i});
        if (length(pathitem) == 2)|| ((length(pathitem) == 3)&&(isempty(pathitem{3})))  %nodes
            [~, st] = ismember(pathitem{1}, NodeTypeStr);
            if st
                j = j+1;
                node(j).id=j;
                GeneIds = checkGeneList(pathitem{2}, GeneDataBase, 0);
                if GeneIds.gid
                    node(j).geneId = GeneIds.gid; node(j).hgncid = GeneIds.hgncid; node(j).ncbiid = GeneIds.ncbiid;
                else
                    node(j).geneId = 0; node(j).hgncid = 0; node(j).ncbiid = 0;
                end
                node(j).name = pathitem{2};
                node(j).type = st;
                node(j).toNode = []; node(j).toType = []; node(j).fromNode = []; node(j).fromType =[];
                fprintf('i:%d      NodeID: [%d, %d]      Type:%d      Name:%s \n', i, j, node(j).geneId, node(j).type, node(j).name);
                
                if st == 2 || st == 4  %number of gene family or complex members
                    node(j).nMembers = 0;
                end
                
            else
                fprintf('Invalid Node Type: [%s %s]\n', pathitem{1},pathitem{2}); pause; pause;
            end
        end
    end
    nNodes=j;

    %*************************************%
    %   Find Edges                        %
    %*************************************%
    fprintf('\n\n Importing Edges\n');
    j = 0;
    for i = 1: nItems
        pathitem = strsplit(TXT{i});
        if (length(pathitem) == 3) && (~isempty(pathitem{3}))
            [~,t] = ismember(pathitem{3}, EdgeTypeStr);
            s=0; d=0;
            for k=1:nNodes
                if strcmp(pathitem{1}, node(k).name)
                    s = k; break;
                end
            end
            for k=1:nNodes
                if strcmp(pathitem{2}, node(k).name)
                    d = k; break;
                end        
            end

            if (s && d && t)  %valid source dest and connection
                j = j+1;
                edge{j} = [s,d,t];

                if mod(i,1000)==0
                    fprintf('i:%d Str:[%s  %s  %s]  Params:[%d - %d - %d] \n',i,  pathitem{1}, pathitem{2}, pathitem{3}, s,d,t);
                end

                node(s).toNode = [node(s).toNode, d];
                node(s).toType = [node(s).toType, t];

                node(d).fromNode = [node(d).fromNode, s];
                node(d).fromType = [node(d).fromType, t];
                
                if ((t == 6) || (t == 5)) && ((node(d).type == 2) || (node(d).type == 4))  %number of gene family or complex members
                    node(d).nMembers = node(d).nMembers + 1;
                end

            else
                fprintf('Error i:%d Str:[%s  %s  %s]  Params:[%d - %d - %d] \n',i,  pathitem{1}, pathitem{2}, pathitem{3}, s,d,t); pause; pause;
            end

        end
        nEdges = j;
    end

    save data3temp
    save(pathway_mat_file, 'TXT',  'node', 'edge', 'nEdges', 'nNodes', 'EdgeTypeStr', 'NodeTypeStr');
    save(pathway_mat_file,'GeneDataBase', '-append'); %for refernce
end



%%
if ipmode == generate_adjc_matrix  %generate agency matrix
    load(pathway_mat_file,  'node');%, 'TXT', 'edge', 'nEdges', 'nNodes', 'EdgeTypeStr', 'NodeTypeStr'); 
    %load('tempnode', 'node');
   
    ProtIds = [];
    for i = 1 : length(node)
        if node(i).type == 3
            ProtIds = [ProtIds, i];
        end
    end
    fprintf('Total Number of Prot Nodes : %d \n', length(ProtIds));
    
    fp=fopen('log_pathway.txt', 'w+'); t=clock; fprintf(fp, 'Date:%d-%d-%d Time:%d:%d:%d \n', t(1:6)); 
    
%     bridge_node_allscenario= {      [1,2,4,5,6],        [2,4,5,6],          [2,4]};
%     VUNodeList_allscenario= {'VUNodeList_12456','VUNodeList_2456','VUNodeList_24'};
%     PNodeList_allscenario = { 'PNodeList_12456', 'PNodeList_2456', 'PNodeList_24'};
% 
%     VDNodeList_allscenario= {'VDNodeList_12456','VDNodeList_2456','VDNodeList_24'};
%     CNodeList_allscenario = { 'CNodeList_12456', 'CNodeList_2456', 'CNodeList_24'};  %children gene list
%     CNodeList2_allscenario = { 'CNodeList2_12456', 'CNodeList2_2456', 'CNodeList2_24'};  %children gene list

    
    bridge_node_allscenario= {[2,3, 4,5,6],   [2,3, 4,5,6],   [2,3, 4,5,6],   [2,3, 4,5,6]};
    MaxLev = [30,   10,  5, 4];
    if learn_upstream
        VUNodeList_allscenario = {'VUNodeList_OTF3_23456_30',  'VUNodeList_OTF3_23456_10',  'VUNodeList_OTF3_23456_5',  'VUNodeList_OTF3_23456_4'};    %includes a correction with respect to VUNodeList_23456
        PNodeList_allscenario  = {'PNodeList_OTF3_23456_30',   'PNodeList_OTF3_23456_10',   'PNodeList_OTF3_23456_5',   'PNodeList_OTF3_23456_4'};
        PNodeList2_allscenario = {'PNodeList_OTF4_23456_30',   'PNodeList_OTF4_23456_10',   'PNodeList_OTF4_23456_5',   'PNodeList_OTF4_23456_4'};  %remove far genes
        
        CNodeListByParrentList  = {'CbyPNodeList_OTF3_23456_30',  'CbyPNodeList_OTF3_23456_10',  'CbyPNodeList_OTF3_23456_5', 'CbyPNodeList_OTF3_23456_4'}
        CNodeListByParrentList2 = {'CbyPNodeList_OTF4_23456_30',  'CbyPNodeList_OTF4_23456_10',  'CbyPNodeList_OTF4_23456_5', 'CbyPNodeList_OTF4_23456_4'}
    end
    if learn_downstream
        %error('not developed for the new version');
        VDNodeList_allscenario= {'VDNodeList_OTF3_23456_30',  'VDNodeList_OTF3_23456_10', 'CNodeList_OTF3_23456_5', 'CNodeList_OTF3_23456_4'};
        CNodeList_allscenario = { 'CNodeList_OTF3_23456_30',  'CNodeList_OTF3_23456_10',  'CNodeList_OTF3_23456_5', 'CNodeList_OTF3_23456_4'};  %children gene list
        CNodeList2_allscenario = {'CNodeList_OTF4_23456_30',  'CNodeList_OTF4_23456_10',  'CNodeList_OTF4_23456_5', 'CNodeList_OTF4_23456_4'};  %children gene list, based on parent list
    end    
    
    
    %make parent list and depth for each protein
    for scenario = [3,2]
        browse_params.debug_active=0;
        browse_params.bridge_node = bridge_node_allscenario{scenario};
        browse_params.term_node = 3;
        browse_params.MaxLev = MaxLev(scenario);
        browse_params.node = node;
        
        if learn_upstream
            browse_params.valid_edge_type_first_ring        = [3,4];         %trnscript:[3:4]
            browse_params.valid_edge_type_intermediate_gene = [1,2,   5,6];
            
            VUNodeList = VUNodeList_allscenario{scenario};
            PNodeList  = PNodeList_allscenario{scenario};
            PNodeList2  = PNodeList2_allscenario{scenario};
            CbyPNodeList   = CNodeListByParrentList{scenario};
            CbyPNodeList2  = CNodeListByParrentList2{scenario};
            childList2={};
        end        
        if learn_downstream
            browse_params.valid_edge_type_term_node = [3,4];
            
            VDNodeList = VDNodeList_allscenario{scenario};
            CNodeList  = CNodeList_allscenario{scenario};
            CNodeList2 = CNodeList2_allscenario{scenario};
        end
        
        
        fprintf('\nScenario :%d/%d \n', scenario, length(MaxLev)); gids=[];
        for i=1:length(ProtIds)
            if ~mod(i,100), fprintf('.'); end
            nodeid = ProtIds(i); rootid = ProtIds(i); 
            gid = node(nodeid).geneId;  
            if gid > 0
                %sweep Upstream
                if learn_upstream
                    visitUpList = []; parentList=[];  smatrix = []; depth=0;  ConType=0; 
                    %[~, visitUpList, parentList] = visit_nodeU(nodeid, rootid, depth,  term_node, bridge_node, node, , MaxLev(scenario), 0);
                    [~, visitUpList, parentList] =  visit_nodeU(nodeid, rootid, depth,  visitUpList, parentList,smatrix, ConType, browse_params);
                    
                    
                    debug_active = browse_params.debug_active;
% term_node = browse_params.term_node;
% bridge_node = browse_params.bridge_node;
% valid_edge_type_term_node = browse_params.valid_edge_type_term_node;
% node= browse_params.node;
% MaxLev = browse_params.MaxLev;
%valid_edge_type_term_node = [3,4];   %set to transcription:[3,4]
                    
                    for j = 1:size(parentList,2)
                        parentList(4,j) = node(parentList(1,j)).geneId;  %row1:node id,   row2:depth,   row3:first link connection Type    row4: geneId
                    end
                    if ~isempty(parentList), gids=[gids,gid]; end
                    parentList = remove_double_parents(parentList,0);
                    parentList2 = remove_double_parents(parentList,1);
                    eval([VUNodeList , '{gid} = visitUpList;']);
                    eval([PNodeList , '{gid} = parentList;']);
                    eval([PNodeList2 , '{gid} = parentList2;']);
                end
                
                if learn_downstream
                    %sweep DownStream
                    visitDownList = []; childList=[];  smatrix = []; depth=0; ConType=0; UpDir = 0;
                    %[depth, visitDownList, childList] = visit_node(UpDir, nodeid, rootid, depth,  ConType, term_node, bridge_node, node, visitDownList, childList, smatrix, MaxLev(scenario), 0);
                    [~, visitDownList, childList0] = visit_nodeD(nodeid, rootid, depth,  visitDownList, childList, smatrix, ConType, browse_params);

                    
                    for j = 1:size(childList0,2)
                        childList0(4,j) = node(childList0(1,j)).geneId;   %row1:node id,   row2:depth,   row3:first link connection Type    row4: geneId
                    end
                    childList = remove_double_parents(childList0, 0);  %remove double genes
                    childList2 = remove_double_parents(childList0, 1); %remove further genes
                    if sum(sum(childList-childList0)) ~=0 , fprintf('Warning repeated gene in the list: %d %d \n', size(childList0,2), size(childList,2));end
                    
                    eval([VDNodeList , '{gid} = visitDownList;']);
                    eval([CNodeList , '{gid} = childList;']);
                    eval([CNodeList2 , '{gid} = childList2;']);
                end
                

                %debug
                if 0
                    fprintf(fp, '\n******************************************************\n');
                    fprintf(fp, '[%d/%d]:  Viditing node:%d  geneid:%d   gene:%s \n', i, length(ProtIds), nodeid, gid, node(nodeid).name );
                    fprintf(fp, 'Visiting List:%s \n', num2str(visitList(1,:)));
                    for j=1:length(visitList)
                        fprintf(fp, '[%s]  ', node(visitList(1,j)).name);
                    end

                    fprintf(fp, '\nParents List:%s\n', num2str(parentList(1,:)));
                    for j=1:length(parentList)
                        fprintf(fp, '[%s]  ', node(parentList(1,j)).name);
                    end
                end

%                 if learn_upstream && ~isempty(parentList)
%                     for j=1:size(parentList,2)
%                         parentNode = parentList(1,j); depth =  parentList(2,j); parentGene =  parentList(3,j);
%                         if (parentGene>0) && (parentGene <= length(childList2))
%                             childList2{parentGene} = [childList2{parentGene}, [nodeid; depth; node(i).geneId]];
%                         elseif (parentGene>0) && (parentGene > length(childList2))
%                             childList2{parentGene} = [nodeid; depth; node(i).geneId];
%                         end
%                     end
%                 end
            end
        end 

        
        
        
        if learn_upstream 
            fprintf('Building Child List from Parent List ....\n');
            eval(['PN = ', PNodeList ,';']);   
            CN(1:22000)={[]};
            for i = 1: length(PN)
                if size(PN{i},2)> 0  %not empty
                    child = i; child_nodeid = 0;  for k=1:length(node), if node(k).geneId== child, child_nodeid = k; end; end
                    for j=1:size(PN{i},2)
                        parent = PN{i}(4,j);
                        if parent > 0 && parent < 24000 
                            if parent > length(CN), CN{parent}=[]; end
                            if isempty(CN{parent})
                                CN{parent} = [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child];
                            else
                                CN{parent} = [CN{parent}, [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child]];
                            end
                        end
                    end
                end
            end
            eval([CNodeListByParrentList , ' = CN;']);
            
            
            eval(['PN = ', PNodeList2 ,';']);   
            CN(1:22000)={[]};
            for i = 1: length(PN)
                if size(PN{i},2)> 0  %not empty
                    child = i; child_nodeid = 0;  for k=1:length(node), if node(k).geneId== child, child_nodeid = k; end; end
                    for j=1:size(PN{i},2)
                        parent = PN{i}(4,j);
                        if parent > 0 && parent < 24000 
                            if parent > length(CN), CN{parent}=[]; end
                            if isempty(CN{parent})
                                CN{parent} = [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child];
                            else
                                CN{parent} = [CN{parent}, [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child]];
                            end
                        end
                    end
                end
            end
            eval([CNodeListByParrentList2 , ' = CN;']);
            fprintf('Building Child List from Parent List (Finished).\n');
        end
        %Later form child list from parent list to verify

        
        %CHECK WHY CNodeList2 ~= CNodeList 
        
        if learn_upstream
            eval(['save(pathway_mat_file, ''ProtIds'', ''MaxLev'', ''', VUNodeList, ''',  ''', PNodeList, ''', ''', PNodeList2, ''', ''', CNodeList2, ''', ''-append'');']);
        end
        if learn_downstream
            eval(['save(pathway_mat_file, ''ProtIds'', ''MaxLev'', ''', VDNodeList, ''', ''', CNodeList, ''', ''-append'');']);
        end
        
    end
    fclose(fp);
end
    
%%
if ipmode == show_data  %generate agency matrix
    load(pathway_mat_file, 'GeneDataBase', 'TXT',  'node', 'edge', 'nEdges', 'nNodes', 'EdgeTypeStr', 'NodeTypeStr', 'AdjMatrix', 'VUNodeList_12456','PNodeList_12456', 'PNodeList_23456', 'VDNodeList_12456','CNodeList_12456', 'PNodeList_OTF3_23456_5' );    
    GeneList = GeneDataBase.GeneList; ng= length(GeneList); nEdgeType=length(EdgeTypeStr); nNodeType = length(NodeTypeStr);
      
    show_all_nodes = 0;
    show_all_edges_types = 0;
    show_hotspot_nodes = 0;
    show_parent_numbers=1;
    
    if show_all_nodes
        for nt = 1:6
            for i=1:length(node)
                if node(i).type == nt
                   fprintf('i:%d      NodeID: [%d]      Type:%d      Name:%s \n', i, node(i).geneId, node(i).type, node(i).name);
                end
            end
            fprintf('\n***********************************\n\n'); pause;
        end
    end
    
    
    if show_all_edges_types
        EdgeTypeStrNoDir = {'-a', '-a', '-t', '-t', 'component>', 'member>'}; 
        AllConnectionType = {}; AllConnectionCnt = []; ii=0; C = [];
        for i=1: length(edge)
            if isnumeric(edge{i}(1)) 
                sid = edge{i}(1); did = edge{i}(2); ctype = edge{i}(3); stype=node(sid).type; dtype=node(did).type; 
                ConnectionType = sprintf('%s==(%s)==>%s', NodeTypeStr{stype}, EdgeTypeStrNoDir{ctype}, NodeTypeStr{dtype});
                [a,b] = ismember(ConnectionType, AllConnectionType);
                if a,  %existing type
                    AllConnectionCnt(b)=AllConnectionCnt(b)+1;
                else   %new type
                    ii=ii+1; AllConnectionType{ii} = ConnectionType; AllConnectionCnt(ii)=1; C(ii).stype= stype; C(ii).dtype= dtype; C(ii).ctype= ctype; 
                end
            end
        end
        
        for i=1:length(AllConnectionCnt)
              fprintf('%s  \t  Count:%d\n', AllConnectionType{i}, AllConnectionCnt(i));
        end
        fprintf('\n***********************************\n\n\n\n'); pause;
        
        %categorize per source type
        for stype = 1: 6
            for i=1:length(AllConnectionCnt)
                if C(i).stype == stype
                    fprintf('%s  \t  Count:%d\n', AllConnectionType{i}, AllConnectionCnt(i));
                end
            end
            fprintf('\n***********************************\n\n'); pause;
        end        
        
        fprintf('\n\n\n\n');
        %categorize per dest type
        for dtype = [3]%1: 6
            for i=1:length(AllConnectionCnt)
                if C(i).dtype == dtype
                    fprintf('%s  \t  Count:%d\n', AllConnectionType{i}, AllConnectionCnt(i));
                end
            end
            fprintf('\n***********************************\n\n'); pause;
        end        
    end

    if show_hotspot_nodes
        for i=1: length(PNodeList_12456)
            len1(i) = length(PNodeList_12456{i});
        end
        for i=1: length(CNodeList_12456)
            len2(i) = length(CNodeList_12456{i});
        end
        
        ind1= find(len1>0); ind2= find(len2>0);         ind=intersect(ind1,ind2);
        ind1h = find(len1>5); ind2h = find(len2>5);     indh=intersect(ind1h,ind2h);
        ind1hh = find(len1>25); ind2hh = find(len2>25); indhh=intersect(ind1hh,ind2hh);
        
        for i=1:length(indhh)
            fprintf(' g[%d]:%s  ', indhh(i), GeneList{i});
        end
    end
    
    if show_parent_numbers
        for OTF = [0,1,2, 4]
            n1=0; nt=0;
            if OTF==4
                PNodeList = PNodeList_OTF3_23456_5;
            elseif OTF==2
                PNodeList = PNodeList_23456;
            else
                PNodeList = PNodeList_12456;
            end

            ng=min(length(GeneList), length(PNodeList));
            vecLen = zeros(1,ng);
            for gid =1 : ng
                P = PNodeList{gid}; 
                if ~isempty(P)
                        P=P(:,P(4,:)>0);
                    if OTF > 0
                        P=P(:,P(3,:)==3 | P(3,:)==4);
                    end
                    vecLen(gid)=size(P,2);
                end
                if vecLen(gid)> 700
                    %fprintf('gid:%d  gene:%s   Number of parents:%d \n',  gid, GeneList{gid}, vecLen(gid));
                end
                if vecLen(gid) == 1, n1=n1+1; end
                if vecLen(gid) >0, nt=nt+1; end
                %if ~mod(gid,50), pause; end
            end
            fprintf('OTF:%d  %d genes have only one parent\n',OTF, n1);
            fprintf('OTF:%d  %d genes have parents\n', OTF, nt);
            [~,gid] = max(vecLen); 
            fprintf('OTF:%d  gid:%d  gene:%s   Number of parents:%d \n', OTF, gid, GeneList{gid}, vecLen(gid))
            
            x = sort(unique(vecLen)); h = histc(vecLen, x); bar(x,h); title('Gene Parent Count Histogram');% xlim([-10 800]); ylim([-10 120]);
            xlabel('Number of ancestor genes'); ylabel('Histogram (in log scale)'); title(''); saveas(gcf, ['HsitGene_OTF', num2str(OTF), '.fig']); pause;
        end
    end
    
end    
    

%%
if ipmode == browse_node || ipmode == show_network
    
    DirUpv = [0, 1];
    OTF = 4;   
    %OTF 0 browse all links, stop if reaches a gene
    %OTF 1 browse only trnscriptional links, stop if reaches a gene
    %OTF 2 browse only trnscriptional links, pass genes, stop if reaches max depth
    %OTF 3 browse only trnscriptional links, pass genes, stop if reaches max depth
    %OTF 4 =3, limit the number of parents/childs
    
    %UPSTRTEAM PARMAS FOR OTF=2,3 and 4 
    browse_params.bridge_node=[2:6];% Important, bridge nodes include gene
    browse_params.valid_edge_type_first_ring        = [3,4]% only for test later back to [3,4,5,6];         %trnscript:[3:4]
    browse_params.valid_edge_type_intermediate_gene =  [2, 5,6]%only for test later back to [2, 3,4  5,6];  
    browse_params.term_node = 3;
    
    browse_params.valid_edge_type_term_node = [3,4];% [3,4];  %downlink
    
    
    browse_params.MaxLev = 4;
    browse_params.debug_active = 0;
    
    %term_node = browse_params.term_node;
    %MaxLev = browse_params.MaxLev;
    debug_active = 1;
    
    
    if ~exist('edge', 'var'), load(pathway_mat_file,  'node', 'edge', 'EdgeTypeStr', 'NodeTypeStr', 'ProtIds'); end
    if ~exist('GeneDataBase', 'var'), load (GeneNameListFile, 'GeneDataBase'); end
    UpDownStr = {'DownStream', 'UpStream'}; 
   
    browse_params.node = node;
    nodeid = []; gene = input('Enter gene Name or gene Id  [-] or NCBI Id: ', 's'); 
    if isempty(gene)
        
        %from tumorportal: 
        
        gHighlyMutated_BCRA = {'PIK3CA', 'TP53', 'GATA3', 'MAP3K1', 'MLL3', 'CDH1', 'NCOR1', 'MAP2K4', 'PTEN', 'RUNX1', 'PIK3R1', 'CTCF', 'AKT1', 'CBFB', 'SPEN', 'SF3B1', 'ARID1A', 'RB1', 'MLL', 'KRAS'};
        gSignificantlyMutated_BCRA = {'TBX3', 'ERBB2', 'FOXA1', 'MED23', 'STAG2', 'MYB', 'TBL1XR1', 'HIST1H3B', 'CASP8', 'CDKN1B', 'CUL4B', 'RAB40A'};
        gLessSignificantlyMutated_BCRA = {'ERBB3', 'CDC42BPA', 'SETDB1', 'FGFR2', 'GNPTAB', 'EP300', 'ACVR1B'};
        glist = {'PIK3CA', 'TP53', 'GATA3', 'MAP3K1', ...
            'CTNNB1',...
            'MLL3', 'CDH1', 'NCOR1', 'MAP2K4', 'PTEN', 'RUNX1', 'PIK3R1', 'CTCF', 'AKT1', 'CBFB', 'SPEN', 'SF3B1', 'ARID1A', 'RB1', 'MLL', 'KRAS', ... %HIGHLY MUTATED
            'TBX3', 'ERBB2', 'FOXA1', 'MED23', 'STAG2', 'MYB', 'TBL1XR1', 'HIST1H3B', 'CASP8', 'CDKN1B', 'CUL4B', 'RAB40A', ...   %SIGNIFICANTLY MUTATED
            'ERBB3', 'CDC42BPA', 'SETDB1', 'FGFR2', 'GNPTAB', 'EP300', 'ACVR1B'};
        glist = {'APC', 'SMO', 'AR', 'CASP3','APC', 'PSA','APC', 'DCC', 'PKC', 'PKA', 'PI3K', 'DAG', 'HPH'}; 
        glist = GeneDataBase.GeneList(15000:20000)
        glist = GeneDataBase.GeneList(15263);
        
        pause;
    else
        glist  = {gene};
    end
    
    for gidinlist = 1:length(glist)
        gene = glist{gidinlist}; nodeid=[];
        if ~isempty(gene)
            if ~isempty(str2num(gene)), gidval=str2num(gene);
                if str2num(gene) > 0  %NCBI ID
                    [GeneIds] = return_gene_id(gidval, 'ncbi', GeneDataBase);
                else
                    [GeneIds] = return_gene_id(-gidval, 'gid', GeneDataBase);
                end
            else
                [GeneIds] = return_gene_id(gene, 'name', GeneDataBase);
                %[gid, ncbiid, hgncid] = return_gene_id(gene, 'name', GeneDataBase);
            end    
            gid=GeneIds.gid; 
            if gid
                for i=1: length(node)
                    if node(i).type == 3 && node(i).geneId==gid, nodeid=i; end
                end
                fprintf('\nNodeId:%d   GeneId:%d  NCBI:%d   HGNC:%d  gene:%s \n', nodeid, gid,  GeneIds.ncbiid, GeneIds.hgncid, GeneIds.name);
            end
        end

        if isempty(nodeid) || ~gid
            fprintf('No information available for this gene:%d [%s]\n\n', gid, gene);
        else
            rootid=nodeid; 
            ParentChildStr = {'Children', 'Parents'};
            
            for DirUp = DirUpv
                
                if ipmode == browse_node && (ismember(0, OTF) || ismember(1, OTF)) %browse until reach one gene [OTF 0 browse all links, 1: only transcriptional nodes');
                    browse_params.bridge_node = [2,4,5,6];
                    
                    visitList = []; parentList=[]; smatrix=[]; depth = 0; ConType = 0; debug_active=0;
                    [depth , visitList, parentList, smatrix, ConType] = visit_node(DirUp, nodeid, rootid, depth, ConType, browse_params.term_node, bridge_node, node, visitList, parentList, smatrix, browse_params.MaxLev, debug_active);

                    if show_all_parents
                        fprintf('\n******************************************************\n');
                        fprintf('[%d/%d]:  Viditing node:%d  geneid:%d   gene:%s \n', i, length(ProtIds), nodeid, gid, node(nodeid).name );
                        fprintf('\nVisiting [%s] List:%s \n', UpDownStr{DirUp+1}, num2str(visitList)); 
                        for j=1:size(visitList,2), fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, visitList(1,j), visitList(2,j), node(visitList(1,j)).name); end; pause;

                        fprintf('\n\n\n Termination [%s] List:%s\n', ParentChildStr{DirUp+1}, num2str(parentList));
                        for j=1:size(parentList,2), fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, parentList(1,j), parentList(2,j), node(parentList(1,j)).name); end; pause;
                        fprintf('Total Parenets for g[%d]:%s is %d \n', gid, gene, length(parentList));
                    end
                end
                

                if ipmode == browse_node && (ismember(2, OTF) || ismember(3, OTF) || ismember(4, OTF))   %new version , pass after reaching a gene           
                    
                    
                    visitList2 = []; parentList2=[]; smatrix2=[]; depth = 0; ConType = 0; debug_active=0;
                    if DirUp == 1
                            [depth2 , visitList2, parentList2, smatrix2, ConType2] = visit_nodeU(nodeid, rootid, depth,  visitList2, parentList2, smatrix2, ConType, browse_params);

                    else
                            [depth2 , visitList2, parentList2, smatrix2, ConType2] = visit_nodeD(nodeid, rootid, depth, visitList2, parentList2, smatrix2, ConType, browse_params);
                    end
                    
                    
                    if show_all_parents
                        fprintf('\n******************************************************\n');
                        fprintf('[%d/%d]:  Viditing node:%d  geneid:%d   gene:%s \n', i, length(ProtIds), nodeid, gid, node(nodeid).name );
                        fprintf('\nVisiting [%s] List:%s \n', UpDownStr{DirUp+1}, num2str(visitList2)); 
                        for j=1:size(visitList2,2), fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, visitList2(1,j), visitList2(2,j), node(visitList2(1,j)).name); end; pause;

                        fprintf('\n\n\n Termination [%s] List:%s\n', ParentChildStr{DirUp+1}, num2str(parentList));
                        for j=1:size(parentList2,2), fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, parentList2(1,j), parentList2(2,j), node(parentList2(1,j)).name); end; pause;
                        fprintf('Total Parenets for g[%d]:%s is %d \n', gid, gene, length(parentList2));
                    end
                    
                end
                
                
                %PLOT for Browse  node
                %fprintf('Number of %s for Gene [%d]: %s is %d\n', ParentChildStr{DirUp+1}, gid, node(nodeid).name, length(parentList));
                
                if ipmode == browse_node && exist('smatrix', 'var')
                    if isempty(smatrix)
                        fprintf('OTF1: No %s for Gene [%d]: %s \n', ParentChildStr{DirUp+1}, gid, node(nodeid).name);
                    else
                        mygview(DirUp, smatrix, node, 400);
                        [handlelast, handle] = getg(); set(handlelast, 'name', ['Graph[', UpDownStr{1+DirUp},'] for Gene: ', gene]); 
                    end
                end
                if ipmode == browse_node && exist('smatrix2', 'var')
                    if isempty(smatrix2)
                        fprintf('OTF2: No %s for Gene [%d]: %s \n', ParentChildStr{DirUp+1}, gid, node(nodeid).name);
                    else

                        %SAMPLE PLOT FOR PAPER
%                         smatrix2 = remove_extra_nodes_for_plot(smatrix2);
%                         indices_for_OTF4_g15263 = [1:17, 22, 41, 42, 43,44,49,50];  %sparse representation for gene:15263
%                         if gid == 15263,
%                             smatrix2 = smatrix2(indices_for_OTF4_g15263,:); 
%                         end
%                             
                        mygview(DirUp, smatrix2, node, 400);
                        [handlelast, handle] = getg(); set(handlelast, 'name', ['Graph[', UpDownStr{1+DirUp},'] for Gene: ', gene]); 
                    end
                end
            end
            
            
            
            
            
            
            
            %%
           %SHOW PARENTS AND CHILDREN 
           if ipmode == show_network
                if ~exist('PNodeList_OTF3_23456_10','var') || ~exist('PNodeList','var') || ~exist('CNodeList','var')  
                    load(pathway_mat_file, 'PNodeList_OTF3_23456_10', 'PNodeList_OTF4_23456_10', 'VUNodeList_OTF3_23456_10', 'VDNodeList_OTF3_23456_10', 'CbyPNodeList_OTF3_23456_10', 'CNodeList_OTF3_23456_10', 'node', 'edge', 'nEdges', 'nNodes');
                    PNodeList = PNodeList_OTF4_23456_10;
                    VUNodeList = VUNodeList_OTF3_23456_10;
                    VDNodeList = VDNodeList_OTF3_23456_10;
                    CNodeList = CbyPNodeList_OTF3_23456_10; %CNodeList_OTF3_23456_10; [Perhaps nore coherent]
                end
                if ~exist('GeneDataBase', 'var'), load (GeneNameListFile, 'GeneDataBase'); end
                L= input('Enter Maximum Depth [Def:5]:'); if isempty(L),L=5; end
                fprintf('\n******************************************************\n');
                fprintf('Geneid:%d   gene:%s \n', gid, GeneDataBase.GeneList{gid});
                fprintf('\n******************************************************\n');

                if 0
                    visitList = VUNodeList{gid};
                    if ~isempty(visitList) && 0
                        fprintf('\nUPstream Visiting List Node Ids [%d]:%s \n', size(visitList,2), num2str(visitList(1,:))); pause;
                        for j=1:size(visitList,2), 
                            if visitList(2,j) <= L, fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, visitList(j), visitList(2,j), node(visitList(1,j)).name); end; 
                        end; pause;
                    end
                end

                parentList = PNodeList{gid};  
                if ~isempty(parentList) && 1
                    [~,ind] = sort(parentList(2,:)); parentList =parentList(:, ind); 
                    fprintf('Show Parents Up to Level:%d for gene:%s \n', L, GeneDataBase.GeneList{gid}); pause; pause;
                    for i = 1: L, ind = (parentList(2,:)==i);  p = numel(unique(parentList(4,ind))); fprintf('Level:%d   Num[%d]     \n', i, p); end
                    fprintf('\n\n\n Parent List GeneIDS[%d]:%s\n', size(parentList,2), num2str(parentList(4,:))); pause;
                    for j=1:size(parentList,2), 
                        if parentList(2,j) <= L,  fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, parentList(1,j), parentList(2,j), node(parentList(1,j)).name); end;
                    end; pause;
                end


                if 0
                    visitList = VDNodeList{gid};
                    if ~isempty(visitList) && 0
                        fprintf('\nDownStream Visiting List Node Ids [%d]:%s \n', size(visitList,2), num2str(visitList(1,:))); pause;
                        for j=1:size(visitList,2), 
                            if visitList(2,j) <= L,  fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, visitList(j), visitList(2,j), node(visitList(1,j)).name); end
                        end; pause;
                    end
                end

                childrenList = CNodeList{gid}; 
                if ~isempty(childrenList) && 1
                    [~,ind] = sort(childrenList(2,:)); childrenList =childrenList(:, ind); 
                    fprintf('Show Children Up to Level:%d for gene:%s \n', L, GeneDataBase.GeneList{gid}); pause; pause;
                    for i = 1: L, ind = (childrenList(2,:)==i);  p = numel(unique(childrenList(4,ind))); fprintf('Level:%d   Num[%d]     \n', i, p); end
                    fprintf('\n\n\n Children List GeneIDS[%d]:%s\n', size(childrenList,2), num2str(childrenList(4,:))); pause;

                    for j=1:size(childrenList,2), 
                        if childrenList(2,j) <= L, fprintf('j:%d  nodeid:%d  depth:%d [%s] \n', j, childrenList(1,j), childrenList(2,j), node(childrenList(1,j)).name); end
                    end; pause;
                end
           end
            
        end
    end    
end
    




    













standalone_learn_children_from_parent_list=0;
if standalone_learn_children_from_parent_list  %temp script for missed files, already included above
%   clear
    initvars; 
    clearvars -except pathway_mat_file
    
    PNodeListV = {'PNodeList_OTF3_23456', 'PNodeList_OTF3_23456_10', 'PNodeList_OTF3_23456_30', 'PNodeList_OTF3_23456_4', 'PNodeList_OTF3_23456_5', 'PNodeList_OTF4_23456_10', 'PNodeList_OTF4_23456_30', 'PNodeList_OTF4_23456_4', 'PNodeList_OTF4_23456_5'};
    for v = 1: length(PNodeListV)
        PNodeListStr = PNodeListV{v}
        CNodeListStr = ['CbyPNodeList_', PNodeListStr(11:end)]
    
 
        eval(['load(pathway_mat_file, ''node'', ''', PNodeListStr, ''');']);
        eval(['PN = ', PNodeListStr, ';']);
        CN(1:22500)={[]};
        for i = 1: length(PN)
            if size(PN{i},2)> 0  %not empty
                child = i; child_nodeid = 0;  for k=1:length(node), if node(k).geneId== child, child_nodeid = k; end; end
                for j=1:size(PN{i},2)
                    parent = PN{i}(4,j);
                    if parent > 0 && parent < 24000 
                        if parent > length(CN), CN{parent}=[]; end
                        if isempty(CN{parent})
                            CN{parent} = [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child];
                        else
                            CN{parent} = [CN{parent}, [child_nodeid ;PN{i}(2,j); PN{i}(3,j) ;child]];
                        end
                    end
                end
            end
        end
        eval([CNodeListStr , ' = CN;']);
        eval(['save(pathway_mat_file, ''', CNodeListStr , ''', ''-append''); ']);
%         who
%         pause;

    end
end



























% % % % % if 0 % old version generate_matrix  %generate agency matrix
% % % % %     load(pathway_mat_file, 'TXT',  'node', 'edge', 'nEdges', 'nNodes', 'EdgeTypeStr', 'NodeTypeStr', 'AdjMatrix');    
% % % % %     load (GeneNameListFile, 'GeneList', 'GeneID_in_GEXP','GeneID_in_RNASeq', 'GeneNameSynonimsList');
% % % % %     %load (GeneSynNameListFile, 'GeneNameSynonimsList');  %call once out of function for fast run
% % % % %     ng= length(GeneList); nEdgeType=length(EdgeTypeStr); nNodeType = length(NodeTypeStr);
% % % % %     
% % % % %     if 0%~exist('AdjMatrix', 'var')
% % % % %         CounterMatrix=zeros(nNodeType);
% % % % % 
% % % % %         for i=1: length(edge)
% % % % %             if isnumeric(edge{i}(1)) 
% % % % %                 sid = edge{i}(1); did = edge{i}(2); ctype = edge{i}(3);
% % % % %                 fprintf('i:%d / %d  Edge:  [%d:%s [Type:%d]        %d:%s [Type:%d]       %d[%s]] ...   \n', i, length(edge), sid, node(sid).name, node(sid).type, did, node(did).name, node(did).type,  ctype, NodeTypeStr{ctype});
% % % % %                 if (sid > 0) && (sid < length(node)) && (did > 0) && (did < length(node))
% % % % %                     %if node(sid).type < size(CounterMatrix,1) &&  node(did).type < size(CounterMatrix,2)
% % % % %                         CounterMatrix(node(sid).type, node(did).type) = 1 + CounterMatrix(node(sid).type, node(did).type);
% % % % %                     %else
% % % % %                     %    CounterMatrix(node(sid).type, node(did).type) = 1;
% % % % %                     %end
% % % % %                 end
% % % % %             end
% % % % %         end
% % % % %         save(pathway_mat_file, 'CounterMatrix', '-append');    
% % % % %     end
% % % % %     
% % % % %     
% % % % %     if 0%~exist('AdjMatrix', 'var')   %old version
% % % % %         AdjMatrix = zeros(ng);
% % % % %         for i=1: length(edge)
% % % % % 
% % % % %             %check sourse of the edge
% % % % %             sGeneId = [];
% % % % %             fprintf('i:%d / %d  Edge:  [%d   %d   %d] ...   ', i, length(edge), edge{i});
% % % % %             if isnumeric(edge{i}(1)) 
% % % % %                 if (edge{i}(1) > 0) && (edge{i}(1) < length(node))
% % % % %                     sGene = node(edge{i}(1)).name;
% % % % %                     if ~isempty(sGene),
% % % % %                         [GeneId, GeneList, GeneID_in_GEXP, GeneID_in_RNASeq, GeneList_updated] = checkGeneList(sGene, GeneList, GeneID_in_GEXP, GeneID_in_RNASeq, GeneNameSynonimsList);
% % % % %                         if (GeneId>0) && (GeneId < ng)  %a valid gene
% % % % %                             sGeneId=GeneId;
% % % % %                         end
% % % % %                     end
% % % % %                 end
% % % % %             end
% % % % % 
% % % % %             %check dest of the edge
% % % % %             dGeneId = [];
% % % % %             if ~isempty(sGeneId) && isnumeric(edge{i}(2)) 
% % % % %                 if (edge{i}(2) > 0) && (edge{i}(2) < length(node))
% % % % %                     dGene = node(edge{i}(2)).name;
% % % % %                     if ~isempty(dGene),
% % % % %                         [GeneId, GeneList, GeneID_in_GEXP, GeneID_in_RNASeq, GeneList_updated] = checkGeneList(dGene, GeneList, GeneID_in_GEXP, GeneID_in_RNASeq, GeneNameSynonimsList);
% % % % %                         if (GeneId>0) && (GeneId < ng)  %a valid gene
% % % % %                             dGeneId=GeneId;
% % % % %                         end
% % % % %                     end
% % % % %                 end
% % % % %             end
% % % % % 
% % % % %             %check if source-dest-relation is a valid triplet
% % % % %             if ~isempty(dGeneId) && ~isempty(sGene) && ~isnan(edge{i}(3)) &&  edge{i}(3) > 0 &&  edge{i}(3)< nEdgeType;
% % % % %                 AdjMatrix(sGeneId, dGeneId) =  edge{i}(3); % add connection
% % % % %                 fprintf('Connection Gene[%d:%s]  to Gene[%d:%s]   Type[%d:%s] is added \n', sGeneId, sGene, dGeneId, dGene, edge{i}(3), EdgeTypeStr{edge{i}(3)});
% % % % %             else
% % % % %                 fprintf('Invalid Entry  Gene[%d:%s]  to Gene[%d:%s]   Type[%d:%s] is skipped !!!\n', sGeneId, sGene, dGeneId, dGene, edge{i}(3), EdgeTypeStr{edge{i}(3)});
% % % % %             end
% % % % %         end
% % % % %         
% % % % %         %save(pathway_mat_file, 'AdjMatrix', 'GeneList', '-append');    
% % % % %     end    
% % % % %     
% % % % % end




