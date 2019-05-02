// COARSENING NODES : double_dummy=weighted_deg

#include "cmpfuncs.h"
//#include "LEDA/basic_graph_alg.h"
#include <time.h>
#include <math.h>

//double A0,V1,E0;

#define COLOR(X)  ((G[X].status==seed) ? 'c' : 'f')


//#define PART_OF_GREATEST_WEIGHT 0.01
double cpu_time;
clock_t prev_clock;


void update_time()
{
  clock_t now = clock();
  cpu_time+=(double)(now - prev_clock)/((double)CLOCKS_PER_SEC)/60.0;
  prev_clock = now;
}

int MAX_LEVEL;

int basis_solutions_passed;

array<int> solutions_on_level;

list<TGraphC> AllLevelsGraphs;
list<edge_array<double> > AllLevelsEdgeCosts;
double fix_edge_relative_weight(TGraphC & G, edge e)
{
  //  return G[e].w;
  if((V_Cycle_Number>0)&&(Params.CURRENT_LEVEL==0))
    {
      node u = G.source(e);
      node v = G.target(e);

      double alpha= (double)V_Cycle_Number/(double)(Params.V_Cycle_Iter_Num-1);
      //            if(V_Cycle_Number==1)
      //              alpha = 1;
      //            else if(V_Cycle_Number==2)
      //              alpha = 2;
     return G[e].w/pow(fabs(G[u].S_value-G[v].S_value), alpha);
    }
  else
    return G[e].w;
}



void find_all_edges(node s, TGraphC & G, TGraphC & H, node_array<double> & edges_accum)
{


  list<node> edges_to;
  node x,y,z;
  edge e,f,g;

  node sH = G[s].ptr_to_coarse;
  //  cerr << "Taking v=" << G[s].initial_id << endl;
  list_item e_it, g_it;
  forall_items(e_it, G[s].iw_links)
    {

      double e_iw = G[s].iw_links[e_it].second();
      x = G[s].iw_links[e_it].first();
      //      cerr << "v's iw link to " << G[x].initial_id << "with iw=" << e_iw << endl;

      forall_adj_edges(f, x)
        {
          double f_w = G[f].w;
          y = second_adj_for_edge(f, x, G);

          //          cerr <<  "Taking y=" << G[y].initial_id << endl;
          //           cerr << "iw link to y" <<  "with w=" << f_w << endl;
          if(y!=s)
            {
              forall_items(g_it, G[y].iw_links)
                {
                  double g_iw = G[y].iw_links[g_it].second();
                  z = G[y].iw_links[g_it].first();

                  if((z!=s)&&(z!=x)&&(G[z].status==seed))
                    {

                      //                                    cerr <<  "Taking" << endl;
                  //                  cerr << "iw link to z" <<  "with iw=" << g_iw << endl;
                      //                      if(((G[s].initial_id==828)&&(G[z].initial_id==337))||
                      //                         ((G[z].initial_id==828)&&(G[s].initial_id==337)))
                      //                                            cerr << G[s].initial_id << "---(" << e_iw << ")---" << G[x].initial_id << "^" << COLOR(x) << "---(" << f_w << ")---" << G[y].initial_id << "^" << COLOR(y) << "---(" << g_iw << ")---" << G[z].initial_id << "^" << COLOR(z) << endl;

                      if(G[s].initial_id<G[z].initial_id)
                        {
                          edges_to.push_back(G[z].ptr_to_coarse);
                          //                      cerr << "da" << endl;
                          edges_accum[G[z].ptr_to_coarse]+=e_iw*f_w*g_iw;
                          //                      cerr << "opsa" << endl;
                        }
                    }
                }
            }
        }
    }

  list_item d_it;
  //  cerr << "adding" << endl;
  //  edges_to.sort();
  //  edges_to.unique();

  forall_items(d_it, edges_to)
    {
      //      cout << H[sH].initial_id << "---" << H[edges_to[d_it]].initial_id << endl;
      //      if(H[sH].initial_id<H[edges_to[d_it]].initial_id)
      //        {
      if(edges_accum[edges_to[d_it]]>0)
        {
          CEdge ne; ne.w = edges_accum[edges_to[d_it]];
          H.new_edge(edges_to[d_it], sH, ne);
          //          cerr << "added " << H[edges_to[d_it]].initial_id << "---" << H[sH].initial_id << " with w=" << ne.w << endl;
          //        }
          edges_accum[edges_to[d_it]] = 0;
        }
    }
  //  cerr << "added\n";
}

void filter_edges(TGraphC & H)
  // double_dummy = adj_edges_w
{
  // Recalculation of filtering thershold
  double filtering_threshold = (double)Params.threshold_edge_weight*pow(0.9,log(Params.E0/(double)H.number_of_edges()));
  if(filtering_threshold > Params.threshold_edge_weight)
    filtering_threshold = Params.threshold_edge_weight;
  if(filtering_threshold < 0.00001)
    filtering_threshold = 0.00001;
  if(Params.CURRENT_LEVEL==0)
    filtering_threshold = Params.threshold_edge_weight;
  cerr << "Recalculated filtering threshold = " << filtering_threshold << endl;
  //////////////////////////////////////////

  edge e,f;
  node u, v;
  forall_edges(e, H)
    {
      u = H.source(e);
      v = H.target(e);
      /*
            if(H[u].double_dummy < H[e].w)
              H[u].double_dummy=H[e].w;
            if(H[v].double_dummy < H[e].w)
              H[v].double_dummy=H[e].w;
      */
      H[u].double_dummy+=H[e].w;
      H[v].double_dummy+=H[e].w;
    }

  long N = H.number_of_nodes();
  if(Params.use_edges_abs_filter == false)
    {
      //      int c=0;
      forall_edges(e, H)
        {
          u = H.source(e);
          v = H.target(e);
          double u_sum = H[u].double_dummy;
          //          double u_avg = u_sum / (double)H.degree(u);
          if(u_sum*filtering_threshold /*Params.rel_filter(N)*/ /*threshold_edge_weight*/>H[e].w)//u_sum/100
            {
              double v_sum = H[v].double_dummy;
              //              double v_avg = v_sum / (double)H.degree(v);
              if(v_sum*filtering_threshold /*Params.rel_filter(N)*//*threshold_edge_weight*/>H[e].w)//v_sum/100
                {
                  //                  if((H[e].w<u_avg/2.0) && (H[e].w<v_avg/2.0))
                    H[e].be_deleted = true;
                  //              cerr << H[e].w << " (" << u_sum << ", " << v_sum << ", " << H.degree(v) << ", " << H.degree(u) << "), ";
                }
            }
        }
      //  cerr << endl;
      //      cerr << "Using relative value filtering : " << c << " edges\n";
    }
  else if(Params.use_edges_abs_filter == true)
    {
      forall_edges(e, H)
        H[e].w = H[e].w/2.0;

      int c=0;
      forall_edges(e, H)
        if(H[e].w<Params.threshold_edge_weight)
          {
            H[e].be_deleted = true;
            c++;
          }
      cerr << "Using absolute value filtering : " << c << " edges\n";
      forall_edges(e, H)
        if(H[e].be_deleted == true)
          H.del_edge(e);
      forall_edges(e, H)
        {
          if(H[e].be_deleted==true)
            {
              cerr << "BE DELETED" << endl;
              exit(1);
            }
        }
      return;
    }

  forall_edges(e, H)
    {
      if(H[e].be_deleted == true)
        {
          H.del_edge(e);
        }
    }
}

void construct_new_graph(TGraphC & G, TGraphC & H)
{
  node v;
  node u;

  //  Params.threshold_coarse_edges = Params.threshold_edge_weight;

  /*
  forall_nodes(v, H)
    {
      u = H[v].G_prev_ptr;
      find_all_edges(u, G, H);


      H[v].double_dummy = 0;
    }
*/
  edge e;
node_array<double> edges_accumulator(H);

  // Adding edges from paths of length 1
/*
  forall_edges(e, G)
    {
      u = G.source(e); v = G.target(e);
      if((G[u].status==seed)&&(G[v].status==seed))
        {
          CEdge h; h.w = 2*G[e].w;
          H.new_edge(G[u].ptr_to_coarse, G[v].ptr_to_coarse, h);
          edges_accumulator[G[u].ptr_to_coarse] = 0;
        }
    }
*/
  forall_nodes(v, H)
    {
      u = H[v].G_prev_ptr;
      find_all_edges/*_from_paths_of_lenght_2_3*/(u, G, H, edges_accumulator);

      H[v].double_dummy = 0;

    }

  //  cerr << "Edges defined" << endl;
  //  cerr << "Removing double edges" << endl;
  //  forall_edges(e, H)
  //    if(H[H.source(e)].initial_id<H[H.target(e)].initial_id)
  //      H.hide_edge(e);

  cerr << "Start filtering with " << H.number_of_edges() << endl;
  clock_t filter_start = clock();
  filter_edges(H);
  cerr << "Only Filtering time : " << (double)(clock()-filter_start)/(double)CLOCKS_PER_SEC << endl;

  //  G.restore_all_edges();

  // edge e;
  //  forall_edges(e, H)
  //    H[e].w = H[e].w/2.0;


  /*
  forall_nodes(v, H)
    {
      H[v].total_edges_weight = 0;
      forall_adj_edges(e, v)
        H[v].total_edges_weight+=H[e].w;
    }
  */

  /*
  double nseg = 50;
  double max=H[H.first_edge()].w;
  double min=H[H.first_edge()].w;
  double sum=0;
  forall_edges(e, H)
    {
      if(H[e].w > max) max = H[e].w;
      if(H[e].w < min) min = H[e].w;
      sum+=H[e].w;
    }
  double avg = sum / (double)H.number_of_edges();
  double seg = (max-min)/nseg;

  leda::string filename = i2s(Params.CURRENT_LEVEL);
  filename = "stat_e"+filename+".txt";
  FILE * inFile = fopen(filename, "w");
  fprintf(inFile, "max\t%f\n", max);
  fprintf(inFile, "min\t%f\n", min);
  fprintf(inFile, "seg\t%f\n", seg);

  int distr[(int)nseg];
  for(int j=0; j<nseg; j++) distr[j]=0;

  forall_edges(e, H)
    {
      distr[(int)(H[e].w/seg)]++;
      if((int)(H[e].w/seg)==0)
        {
          //          cerr << H[e].p1 << ", " << H[e].p2 << ", " << H[e].p3 << ", " << endl;
        }
    }

  for(int i=0; i<nseg; i++)
    fprintf(inFile, "seg%d\t%d\n", i,distr[i]);

  fclose(inFile);
  //  forall_edges(e, H)
  //    H[e].w = H[e].w/2.0;
  */
}


void calculate_vertex_weights(TGraphC & G, TGraphC & H)
{
  node v, w, u;
  edge e;

  forall_nodes(v, H)
    {
      w = H[v].G_prev_ptr;
      if (H[v].initial_id != G[w].initial_id)
        {
          cerr << "jopa prejopa" << endl;
          exit(1);
        }
      if(G[w].status != seed)
        {
          cerr << "Vertex weight : ERROR in G_prev_ptr" << endl;
          exit(1);
        }

      double edges_to_fine_sum = 0;
      forall_adj_edges(e, w)
        {
          u = second_adj_for_edge(e, w, G);
          if(G[u].status!=seed)
            {
              edges_to_fine_sum+=G[e].iw*G[u].w;
            }
        }
      H[v].w = G[w].w + edges_to_fine_sum;
      //G[w].w = H[v].w;
    }
}

void check_node_weights(TGraphC & G)
{
  //return ;
  double aa=0;
  node v;
  forall_nodes(v, G)
    aa+=G[v].w;
  cerr << "NODES WEIGHT = " << aa << endl;
  return ;
}

void basis_solution_iterations(TGraphC & G)
{
  MAX_LEVEL = 0;
      Params.CURRENT_LEVEL = -1;

      double res = coarsening(G, true);
}

void reset_all_vertices(TGraphC & G)
{
  node v;

  forall_nodes(v, G)
    {
      G[v].status = fine;
    }
}
extern leda::string instance_file;

void mwf_graph_print(TGraphC & G, leda::string filename)
{
  exit(1);
  FILE * outFile;
  outFile = fopen("tempedges.dat", "w");

  cout << filename << endl;
  FILE * allFile;
  allFile = fopen(filename, "w");

  //  TMP_CMP_GRAPHC = &G;
  //  G.sort_nodes(&cmp_ArrId);

    edge e;
    int c=0;
  forall_edges(e, G)
    {
      node v = G.source(e);
      node w = G.target(e);

      if(G[v].ArrId>G[w].ArrId)
	{
	  node t=v; v=w; w=t;
	}

      node t = G.succ_node(v);
      while(t!=w)
	{
	  fprintf(outFile, "%d %d 1.0\n", G[v].initial_id, G[t].initial_id);
	  c++;
	  t = G.succ_node(t);
	}
    }

  fclose(outFile);
  cout << "aga" << endl;
  fprintf(allFile, "p ghct %d %d\n", G.number_of_nodes(), c);
  fclose(allFile);
  //  cout << filename << endl;
  system("cat tempedges.dat  >>"+filename);
  //  system("\\rm -rf tempedges.dat");

}

double V_Cycle_Iterations(TGraphC & G)
{
  node v;
	double res = 0;
  for(int i=0; i<Params.V_Cycle_Iter_Num; i++)
    {
      cerr << "Starting " << i << "-th V-Cycle Iteration" << endl;
      Params.amount_of_coarse_nodes = start_threshold_on_nodes;

      AllLevelsGraphs.clear();
      AllLevelsEdgeCosts.clear();

      reset_all_vertices(G);

      V_Cycle_Number = i;
      //      cerr << "aga" << endl;
      basis_solution_iterations(G);

    }

  return res;
}

void calc_oldpw(TGraphC & G) {
	//  edge_array<double> ret(G);
	edge e;
	node v, u;
	/*
	 forall_edges(e, G)
	 if(((G[G.source(e)].color != green)&&(G[G.target(e)].color != green))||
	 ((G[G.source(e)].color == green)&&(G[G.target(e)].color == green)))
	 G[e].pw = fix_edge_relative_weight(G, e);//G[e].w;
	 */
	forall_edges(e, G)
	G[e].pw = 0;

	forall_nodes(v, G) {
		if (G[v].status!=seed) {
			double edges_to_green_sum = 0;
			forall_adj_edges(e, v) {
				u = second_adj_for_edge(e, v, G);
				if (G[u].status==seed)
					edges_to_green_sum+=fix_edge_relative_weight(G, e);
			}
			forall_adj_edges(e, v) {
				u = second_adj_for_edge(e, v, G);
				if (G[u].status==seed)
					G[e].pw = fix_edge_relative_weight(G, e)/edges_to_green_sum;
			}
		}
	}
}


void coarsening_preparations(TGraphC & G, TGraphC & H)
{
	if((Params.alg_dist==wag_and_gs)||(Params.alg_dist==gauss_seidel)||(Params.iw_calc_type==iw_alg))
	{
		//nodes_algdist_relaxation2d(G);
		//algdist_Jacobi_underrelax(G);
	  /*
	  if(Params.sym_nonsym==symmetric)
	    algdist_symJacobi_underrelax(G);
	  else
	    algdist_nonsymJacobi_underrelax(G);
	  */

	  //algdist_minla_relax(G);
	  algdist_symJacobi_underrelax(G);

	}

	double sum_fines_wd = sort_fine_nodes_by_wdegree(G);

	add_heavy_nodes(G, H, sum_fines_wd);

	if (Params.alg_dist==gauss_seidel)
	{
		coarse_nodes_alg_sum(G, H);
		//wag_coarse_nodes(G, H);
		//coarse_nodes_alg_max(G, H);
		//nodes_algdist_compatible_relaxation2d(G);

	}
	else if(Params.alg_dist==wag_and_gs)
		coarse_nodes_alg_and_wag_sum(G, H);
	else if(Params.alg_dist==jacobian)
	{
		cerr << "Jacobi is not supported" << endl;
		exit(1);
	}
	else if(Params.alg_dist==noalgdist)
	{
		//nodes_algdist_relaxation2d(G);
		//coarse_nodes_alg_sum(G, H);
		wag_coarse_nodes(G, H);
	}

	calc_oldpw(G);
	clear_iw_links(G);

	//if((Params.alg_dist==wag_and_gs)||(Params.alg_dist==gauss_seidel))
	//	nodes_algdist_compatible_relaxation2d(G);

	if(Params.iw_calc_type==iwbamg)
	{
		exit(1);
		if(interpol_order(G)>1)
		{
			//wag_calc_iw(G);
			new_bamg_alg_calc_iw_by_rw(G);
			//bamg_calc_iw(G);
		}
		else
			minla_alg_calc_iw_by_rw(G);
			//alg_calc_iw_by_rw_and_w(G);
			//new_bamg_alg_calc_iw_by_rw(G);
	}
	else if(Params.iw_calc_type==iw_alg)
	{
		//nodes_algdist_compatible_relaxation2d(G);
		alg_calc_iw_by_rw(G, H);
		//alg_calc_iw_by_rw_and_w(G);
	}
	else if(Params.iw_calc_type==iw_alg_and_wag)
		alg_calc_iw_by_rw_and_w(G);
	else if(Params.iw_calc_type==iw_wag)
	{
		//Params.alg_dist=gauss_seidel;
		//nodes_algdist_relaxation2d(G);
		//new_bamg_alg_calc_iw_by_rw(G);
		//alg_calc_iw_by_rw(G);
		//Params.alg_dist=noalgdist;
		wag_calc_iw(G);
	}
}

double coarsening( TGraphC & G,
                  bool first_time_down /* for multiple basis solutions */)

{
  clock_t start, end;

  Params.CURRENT_LEVEL++;
  double avg_deg = print_graph_info(G);

  int level = Params.CURRENT_LEVEL;
  
  if(Params.CURRENT_LEVEL==0) {
	  Params.E0=G.number_of_edges();
	  Params.V0=G.number_of_nodes();
	}


  // Construct smaller graph H
	if(G.number_of_nodes() > Params.min_graph_size)
    {
      TGraphC H;
      clock_t start = clock();
		
	  coarsening_preparations(G, H);

          construct_new_graph(G, H);



          start = clock();
          //          cerr << "Vertex weight : Start" << endl;
          calculate_vertex_weights(G, H);

          amg_hierarchy_output(G, H);


      double H_lin_arr = coarsening(H, first_time_down);

    }
    else {
		TGraphC * Dummy;
		amg_hierarchy_output(G, *Dummy);
	}

  return 0;
}
