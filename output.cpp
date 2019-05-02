#include "cmpfuncs.h"
#include "readprint.h"

const leda::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y_%m_%d_%X", &tstruct);
    for(int i=0; i<80; i++) 
      if(buf[i]==':')
        buf[i] = '_';   
    return buf;
}

void output_visMatlab(TGraphC& G, leda::string fname) {

  leda::string time_now = currentDateTime();

  cerr << "Saving Matlab visualization in " << fname << endl;
  node u, v; 
  edge e;
  
  file_ostream o(fname);
  o << "M = sparse(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  o << "M1s = sparse(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;

  leda::string gr_id_name = "id_"+time_now;
  file_ostream o_gr_id(gr_id_name+".dat");
  
  leda::string gr_arr_name = "arr_"+time_now;
  file_ostream o_gr_arr(gr_arr_name+".dat");
  
  forall_edges(e, G) {
  	u = G.source(e); v = G.target(e);
  	o_gr_id << G[u].initial_id << " " << G[v].initial_id << " " << G[e].w << endl;
  	o_gr_arr << G[u].ArrId << " " << G[v].ArrId << " " << G[e].w << endl;
  }
  o << "load " << gr_id_name << ".dat" << endl;
  o << "load " << gr_arr_name << ".dat" << endl;
  
  o << "M = sparse(" << gr_id_name << "(:,1), " << gr_id_name << "(:,2), " << gr_id_name << "(:,3), " << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  o << "M1s = sparse(" << gr_arr_name << "(:,1), " << gr_arr_name << "(:,2), " << gr_arr_name << "(:,3), " << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
    
  o << "M = M + M';" << endl;
  o << "M1s = M1s + M1s';" << endl;
  o << "figure; spy(M); hold" << endl;
  o << "figure; spy(M1s); hold" << endl;
  o.close();  
  o_gr_id.close();  
  o_gr_arr.close();
//G.sort_nodes(&cmp_ArrId);
}

void output_visMatlab_old(TGraphC& G, leda::string fname) {
  cerr << "Saving Matlab visualization in " << fname << endl;
  node u, v; 
  edge e;
  file_ostream o(fname);
  o << "M = sparse(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  o << "M1s = sparse(" << G.number_of_nodes() << ", " << G.number_of_nodes() << ");" << endl;
  forall_edges(e, G) {
  	u = G.source(e); v = G.target(e);
  	o << "M(" << G[u].initial_id << ", " << G[v].initial_id << ") = " << G[e].w << "; " << "M(" << G[v].initial_id << ", " << G[u].initial_id << ") = " << G[e].w << ";" << endl;
  	o << "M1s(" << G[u].ArrId << ", " << G[v].ArrId << ") = " << G[e].w << "; " << "M1s(" << G[v].ArrId << ", " << G[u].ArrId << ") = " << G[e].w << ";" << endl;
  }
  o << "figure; spy(M); hold" << endl;
  o << "figure; spy(M1s); hold" << endl;
  o.close();  
//G.sort_nodes(&cmp_ArrId);
}


void print_G_lemon(TGraphC& G) {
  leda::string f = "tmp_graph.txt";
  file_ostream o(f);
  node v, w;
  edge e;
  
  double max_T = 0;
  double max_Pi = 0;
  forall_nodes(v, G) {
	if(max_T < G[v].T_max)
	  max_T = G[v].T_max;
	if(max_Pi < G[v].w)
	  max_Pi = G[v].w;
  }
  
  o << "@nodes\n" << "id label T Pi" << endl;
  int i=1;
  forall_nodes(v, G) {
	o << i << " \"" << i << "\" " << G[v].T_max/(2*max_T) << " " << G[v].w/max_Pi << endl;
	G[v].lemon_id = i;
	i++;
  }
  o << "@edges\n"  << "    W" << endl;
  forall_edges(e, G) {
	v = G.source(e); w = G.target(e);
	o  << G[v].lemon_id << " " << G[w].lemon_id << " " << G[e].w << endl;
  }
}


void amg_hierarchy_output(TGraphC& G, TGraphC& H) {
  
	TMP_CMP_GRAPHC = &G;
	G.sort_nodes(&cmp_ArrId);

	leda::string f = "level"+i2s(Params.CURRENT_LEVEL)+".dat";
		
	file_ostream o(f);
		
  	node v, w;
  	edge e;
  	int i=1;
  	

	if(&H!=NULL) {	
		forall_nodes(v, H)
		{
			H[v].ArrId = i;
			i++;
		}
	}
  	o << G.number_of_nodes() << " " << G.number_of_edges() << endl;
	
  	forall_nodes(v, G)
  	{
	
  		if((G[v].status==seed)&&(&H!=NULL))
  			o << "c " << H[G[v].ptr_to_coarse].ArrId << " " << G[v].w << " ";
  		else
  			o << "f " << G[v].w << " ";
  		  	
		forall_adj_edges(e, v) {
  			w = second_adj_for_edge(e, v, G);
  			if(G[v].status == seed)
				o << G[w].ArrId << " " << G[e].w << " ";
			else {
				o << G[w].ArrId << " " << G[e].w << " ";
				if(G[w].status == seed)
					o << "iw:" << G[e].iw << " ";
			}
		}
		
		o << endl;
  	}
}

