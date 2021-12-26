/*
 *   libpal - Automated Placement of Labels Library
 *
 *   Copyright (C) 2008 Maxence Laurent, MIS-TIC, HEIG-VD
 *                      University of Applied Sciences, Western Switzerland
 *                      http://www.hes-so.ch
 *
 *   Contact:
 *      maxence.laurent <at> heig-vd <dot> ch
 *    or
 *      eric.taillard <at> heig-vd <dot> ch
 *
 * This file is part of libpal.
 *
 * libpal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libpal is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libpal.  If not, see <http://www.gnu.org/liceprintCostes/>.
 *
 */

#include "pal.h"
#include "palstat.h"
#include "layer.h"
#include "rtree.hpp"
#include "feature.h"
#include "geomfunction.h"
#include "labelposition.h"
#include "problem.h"
#include "util.h"
#include "priorityqueue.h"
#include "internalexception.h"
#include <cfloat>
#include <limits> //for std::numeric_limits<int>::max()
#include <cassert>
#include <QSet>
#include "qgslabelingengine.h"
//+++++++++++++++++++gpl++++++++++++++++++++++++++++
#include "qgslogger.h"
#include "debugger.h"
#include <unordered_map>
#include <unordered_set>
#include "../../app/aclabeltester.h"
//-------------------gpl----------------------------
using namespace pal;
//+++++++++++debug+++++++++++++++++++++
bool gplDebugger = false;
bool gplPrinter = false;
bool testPrinter = false;
QSet<int> QgsFeatureIDS_prev;
int numF_prev;
int numL_prev;
unordered_map<int, unordered_set<int>> candidates_prev;
//---------------debug-----------------
//*+++++++++++++modification++++++++++++++++++++++++++++++
//TODO: if zooming level or ragion change, the init should be reset to true
//bool init = true;
// solution is cached in it after one round
// map <qgsFeatureID, offset>
//qgsFeatureID: unique id of feature
// offset: the offset of the choosen label to the startID of the feature is unique 
//unordered_map<int, int> solution_prev; 
//--------------modification-----------------------------

inline void delete_chain( Chain *chain )
{
  if ( chain )
  {
    delete[] chain->feat;
    delete[] chain->label;
    delete chain;
  }
}

Problem::Problem()
{
  bbox[0] = 0;
  bbox[1] = 0;
  bbox[2] = 0;
  bbox[3] = 0;
  featWrap = nullptr;
  candidates = new RTree<LabelPosition *, double, 2, double>();
  candidates_sol = new RTree<LabelPosition *, double, 2, double>();
  candidates_subsol = nullptr;
  conflictGraph = nullptr;
}

Problem::~Problem()
{
  if ( sol )
  {
    delete[] sol->s;
    delete sol;
  }

  delete[] featWrap;
  delete[] featStartId;
  delete[] featNbLp;

  qDeleteAll( mLabelPositions );
  mLabelPositions.clear();

  delete[] inactiveCost;

  delete candidates;
  delete candidates_sol;

  delete candidates_subsol;
}

typedef struct
{
  int id;
  double inactiveCost;
  double nbOverlap;
} Ft;

inline bool borderSizeInc( void *l, void *r )
{
  return ( reinterpret_cast< SubPart * >( l ) )->borderSize > ( reinterpret_cast< SubPart * >( r ) )->borderSize;
}

void Problem::reduce()
{

  int i;
  int j;
  int k;

  int counter = 0;

  int lpid;

  bool *ok = new bool[nblp];
  bool run = true;

  for ( i = 0; i < nblp; i++ )
    ok[i] = false;


  double amin[2];
  double amax[2];
  LabelPosition *lp2 = nullptr;

  while ( run )
  {
    run = false;
    for ( i = 0; i < nbft; i++ )
    {
      // ok[i] = true;
      for ( j = 0; j < featNbLp[i]; j++ )  // foreach candidate
      {
        if ( !ok[featStartId[i] + j] )
        {
          if ( mLabelPositions.at( featStartId[i] + j )->getNumOverlaps() == 0 ) // if candidate has no overlap
          {
            run = true;
            ok[featStartId[i] + j] = true;
            // 1) remove worse candidates from candidates
            // 2) update nb_overlaps
            counter += featNbLp[i] - j - 1;

            for ( k = j + 1; k < featNbLp[i]; k++ )
            {

              lpid = featStartId[i] + k;
              ok[lpid] = true;
              lp2 = mLabelPositions.at( lpid );

              lp2->getBoundingBox( amin, amax );

              nbOverlap -= lp2->getNumOverlaps();
              candidates->Search( amin, amax, LabelPosition::removeOverlapCallback, reinterpret_cast< void * >( lp2 ) );
              lp2->removeFromIndex( candidates );
            }

            featNbLp[i] = j + 1;
            break;
          }
        }
      }
    }
  }

  this->nblp -= counter;
  delete[] ok;
}

void Problem::init_sol_empty()
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
if(gplDebugger){
  std::cout<< "init_sol_empty"<<endl;
}
//-------------------gpl-----------------------------------
  int i;

  if ( sol )
  {
    delete[] sol->s;
    delete sol;
  }

  sol = new Sol();
  sol->s = new int[nbft];

  for ( i = 0; i < nbft; i++ )
    sol->s[i] = -1;

  sol->cost = nbft;
}




typedef struct
{
  PriorityQueue *list = nullptr;
  LabelPosition *lp = nullptr;
  RTree <LabelPosition *, double, 2, double> *candidates;
} FalpContext;

bool falpCallback2( LabelPosition *lp, void *ctx )
{
  FalpContext *context = reinterpret_cast< FalpContext * >( ctx );
  LabelPosition *lp2 = context->lp;
  PriorityQueue *list = context->list;

  if ( lp->getId() != lp2->getId() && list->isIn( lp->getId() ) && lp->isInConflict( lp2 ) )
  {
    list->decreaseKey( lp->getId() );
  }
  return true;
}


void ignoreLabel( LabelPosition *lp, PriorityQueue *list, RTree <LabelPosition *, double, 2, double> *candidates )
{


  FalpContext *context = new FalpContext();
  context->candidates = nullptr;
  context->list = list;
  double amin[2];
  double amax[2];

  if ( list->isIn( lp->getId() ) )
  {
    list->remove( lp->getId() );

    lp->getBoundingBox( amin, amax );

    context->lp = lp;
    candidates->Search( amin, amax, falpCallback2, context );
  }

  delete context;
}


bool falpCallback1( LabelPosition *lp, void *ctx )
{
  FalpContext *context = reinterpret_cast< FalpContext * >( ctx );
  LabelPosition *lp2 = context->lp;
  PriorityQueue *list = context->list;
  RTree <LabelPosition *, double, 2, double> *candidates = context->candidates;

  if ( lp2->isInConflict( lp ) )
  {
    ignoreLabel( lp, list, candidates );
  }
  return true;
}



/* Better initial solution
 * Step one FALP (Yamamoto, Camara, Lorena 2005)
 */
void Problem::init_sol_falp(test::Performance& performance, bool final,const bool& is_initial,unordered_map<int, int>& my_solution_prev)
{
  auto start = std::chrono::system_clock::now();
  int solutionCount = 0;
//+++++++++++++++++++gpl++++++++++++++++++++++++++++++++++
if(gplDebugger){
  std::cout<< "init_sol_falp"<<endl;
}
//-------------------gpl-----------------------------------
  int i, j;
  int label;
  PriorityQueue *list = nullptr;

  init_sol_empty();

  list = new PriorityQueue( nblp, all_nblp, true );

  double amin[2];
  double amax[2];

  FalpContext *context = new FalpContext();
  context->candidates = candidates;
  context->list = list;

  LabelPosition *lp = nullptr;

  for ( i = 0; i < nbft; i++ )
    for ( j = 0; j < featNbLp[i]; j++ )
    {
      label = featStartId[i] + j;
      try
      {
        list->insert( label, mLabelPositions.at( label )->getNumOverlaps() );
      }
      catch ( pal::InternalException::Full & )
      {
        continue;
      }
    }

  while ( list->getSize() > 0 ) // O (log size)
  {
    if ( pal->isCanceled() )
    {
      delete context;
      delete list;
      return;
    }

    label = list->getBest();   // O (log size)


    lp = mLabelPositions.at( label );

    if ( lp->getId() != label )
    {
      //error
    }

    int probFeatId = lp->getProblemFeatureId();
    sol->s[probFeatId] = label;
    solutionCount++;

    for ( i = featStartId[probFeatId]; i < featStartId[probFeatId] + featNbLp[probFeatId]; i++ )
    {
      ignoreLabel( mLabelPositions.at( i ), list, candidates );
    }


    lp->getBoundingBox( amin, amax );

    context->lp = lp;
    candidates->Search( amin, amax, falpCallback1, reinterpret_cast< void * >( context ) );
    candidates_sol->Insert( amin, amax, lp );
  }

  delete context;
//+++++++++++++++++++gpl++++++++++++++++++++++++++++++++++
if(gplDebugger){
  solution_cost();
  cout<< sol->cost<<endl;
}
performance.name = "init_sol_falp";
  performance.solutionSize = solutionCount;
 auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
 solution_cost();
 performance.solutionWeight = sol->cost;
  bool testini = is_initial;
 performance.remainingLabels = cacheSolution(final,is_initial, my_solution_prev);
  if(testini) assert(performance.remainingLabels == 0 || performance.remainingLabels == -1);
 if(testPrinter){
   performance.print();
 }
//-------------------gpl-----------------------------------



  if ( displayAll )
  {
    int nbOverlap;
    int start_p;
    LabelPosition *retainedLabel = nullptr;
    int p;

    for ( i = 0; i < nbft; i++ ) // forearch hidden feature
    {
      if ( sol->s[i] == -1 )
      {
        nbOverlap = std::numeric_limits<int>::max();
        start_p = featStartId[i];
        for ( p = 0; p < featNbLp[i]; p++ )
        {
          lp = mLabelPositions.at( start_p + p );
          lp->resetNumOverlaps();

          lp->getBoundingBox( amin, amax );


          candidates_sol->Search( amin, amax, LabelPosition::countOverlapCallback, lp );

          if ( lp->getNumOverlaps() < nbOverlap )
          {
            retainedLabel = lp;
            nbOverlap = lp->getNumOverlaps();
          }
        }
        sol->s[i] = retainedLabel->getId();

        retainedLabel->insertIntoIndex( candidates_sol );

      }
    }
  }

  delete list;
}

void Problem::popmusic(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev)
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
auto start = std::chrono::system_clock::now();
int solutionCount = 0;
if(gplDebugger){
   std::cout<<"popmusic"<<endl;
}
//-------------------gpl-----------------------------------
  if ( nbft == 0 )
    return;

  int i;
  int seed;
  bool *ok = new bool[nbft];

  int r = pal->popmusic_r;

  SearchMethod searchMethod = pal->searchMethod;

  candidates_subsol = new RTree<LabelPosition *, double, 2, double>();

  double delta = 0.0;

  int it = 0;

  SubPart *current = nullptr;

  labelPositionCost = new double[all_nblp];
  nbOlap = new int[all_nblp];

  featWrap = new int[nbft];
  memset( featWrap, -1, sizeof( int ) *nbft );

  SubPart **parts = new SubPart*[nbft];
  int *isIn = new int[nbft];

  memset( isIn, 0, sizeof( int ) *nbft );


  for ( i = 0; i < nbft; i++ )
  {
    parts[i] = subPart( r, i, isIn );
    ok[i] = false;
  }
  delete[] isIn;
  Util::sort( reinterpret_cast< void ** >( parts ), nbft, borderSizeInc );
  //sort ((void**)parts, nbft, borderSizeDec);

  init_sol_falp(performance,false,is_initial,my_solution_prev);
  solution_cost();

  int popit = 0;

  seed = 0;
  while ( true )
  {
    it++;
    /* find the next seed not ok */
    for ( i = ( seed + 1 ) % nbft; ok[i] && i != seed; i = ( i + 1 ) % nbft )
      ;

    if ( i == seed && ok[seed] )
    {
      current = nullptr; // everything is OK :-)
      break;
    }
    else
    {
      seed = i;
      current = parts[seed];
    }

    // update sub part solution
    candidates_subsol->RemoveAll();

    for ( i = 0; i < current->subSize; i++ )
    {
      current->sol[i] = sol->s[current->sub[i]];
      if ( current->sol[i] != -1 )
      {
        mLabelPositions.at( current->sol[i] )->insertIntoIndex( candidates_subsol );
      }
    }

    switch ( searchMethod )
    {
      //case branch_and_bound :
      //delta = current->branch_and_bound_search();
      //   break;

      case POPMUSIC_TABU :
        delta = popmusic_tabu( current );
        break;
      case POPMUSIC_TABU_CHAIN :
        delta = popmusic_tabu_chain( current );
        break;
      case POPMUSIC_CHAIN :
        delta = popmusic_chain( current );
        break;
      default:
        delete[] ok;
        delete[] parts;
        return;
    }

    popit++;

    if ( delta > EPSILON )
    {
      /* Update solution */
      for ( i = 0; i < current->borderSize; i++ )
      {
        ok[current->sub[i]] = false;
      }

      for ( i = current->borderSize; i < current->subSize; i++ )
      {

        if ( sol->s[current->sub[i]] != -1 )
        {
          mLabelPositions.at( sol->s[current->sub[i]] )->removeFromIndex( candidates_sol );
        }

        sol->s[current->sub[i]] = current->sol[i];

        if ( current->sol[i] != -1 )
        {
          mLabelPositions.at( current->sol[i] )->insertIntoIndex( candidates_sol );
        }

        ok[current->sub[i]] = false;
      }
    }
    else  // not improved
    {
      ok[seed] = true;
    }
  }

  solution_cost();

  delete[] labelPositionCost;
  delete[] nbOlap;

  for ( i = 0; i < nbft; i++ )
  {
    delete[] parts[i]->sub;
    delete[] parts[i]->sol;
    delete parts[i];
  }
  delete[] parts;

  delete[] ok;
  solution_cost();
  performance.name = "popmusic";
  performance.solutionWeight = sol->cost;
  for(int i = 0; i < nbft; i++){
    if(sol->s[i] != -1) solutionCount++;
  }
  performance.solutionSize = solutionCount;
  auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
 bool testini = is_initial;
 cout<< "INIT IS********* "<< is_initial << endl;
 performance.remainingLabels = cacheSolution(true,is_initial,my_solution_prev);
  if(testini) assert(performance.remainingLabels == 0);
   if(testPrinter){
   performance.print();
 }
}

typedef struct
{
  QLinkedList<int> *queue;
  int *isIn = nullptr;
  LabelPosition *lp = nullptr;
} SubPartContext;

bool subPartCallback( LabelPosition *lp, void *ctx )
{
  SubPartContext *context = reinterpret_cast< SubPartContext * >( ctx );
  int *isIn = context->isIn;
  QLinkedList<int> *queue = context->queue;


  int id = lp->getProblemFeatureId();
  if ( !isIn[id] && lp->isInConflict( context->lp ) )
  {
    queue->append( id );
    isIn[id] = 1;
  }

  return true;
}

/* Select a sub part, expected size of r, from seed */
SubPart *Problem::subPart( int r, int featseed, int *isIn )
{
  QLinkedList<int> *queue = new QLinkedList<int>;
  QLinkedList<int> *ri = new QLinkedList<int>;

  int *sub = nullptr;

  int id;
  int featS;
  int p;
  int i;

  int n = 0;
  int nb = 0;

  double amin[2];
  double amax[2];

  SubPartContext context;
  context.queue = queue;
  context.isIn = isIn;

  queue->append( featseed );
  isIn[featseed] = 1;

  LabelPosition *lp = nullptr;

  while ( ri->size() < r && !queue->isEmpty() )
  {
    id = queue->takeFirst();
    ri->append( id );

    featS = featStartId[id];
    p = featNbLp[id];

    for ( i = featS; i < featS + p; i++ )  // foreach candidat of feature 'id'
    {
      lp = mLabelPositions.at( i );

      lp->getBoundingBox( amin, amax );

      context.lp = lp;
      candidates->Search( amin, amax, subPartCallback, reinterpret_cast< void * >( &context ) );
    }
  }

  nb = queue->size();
  n = ri->size();

  sub = new int[n + nb];

  i = 0;

  while ( !queue->isEmpty() )
  {
    sub[i] = queue->takeFirst();
    isIn[sub[i]] = 0;
    i++;
  }

  while ( !ri->isEmpty() )
  {
    sub[i] = ri->takeFirst();
    isIn[sub[i]] = 0;
    i++;
  }

  delete queue;
  delete ri;

  SubPart *subPart = new SubPart();

  subPart->probSize = n;
  subPart->borderSize = nb;
  subPart->subSize = n + nb;
  subPart->sub = sub;
  subPart->sol = new int [subPart->subSize];
  subPart->seed = featseed;
  return subPart;
}

double Problem::compute_feature_cost( SubPart *part, int feat_id, int label_id, int *nbOverlap )
{
  double cost;
  *nbOverlap = 0;

  LabelPosition::CountContext context;
  context.inactiveCost = inactiveCost;
  context.nbOv = nbOverlap;
  context.cost = &cost;

  double amin[2];
  double amax[2];
  LabelPosition *lp = nullptr;

  cost = 0.0;

  if ( label_id >= 0 ) // is the feature displayed ?
  {
    lp = mLabelPositions.at( label_id );

    lp->getBoundingBox( amin, amax );

    context.lp = lp;
    candidates_subsol->Search( amin, amax, LabelPosition::countFullOverlapCallback, reinterpret_cast< void * >( &context ) );

    cost += lp->cost();
  }
  else
  {
    cost = inactiveCost[part->sub[feat_id]];
    //(*nbOverlap)++;
  }

  return cost;

}

double Problem::compute_subsolution_cost( SubPart *part, int *s, int *nbOverlap )
{
  int i;
  double cost = 0.0;
  int nbO = 0;

  *nbOverlap = 0;

  for ( i = 0; i < part->subSize; i++ )
  {
    cost += compute_feature_cost( part, i, s[i], &nbO );
    *nbOverlap += nbO;
  }

  return cost;
}



typedef struct _Triple
{
  int feat_id;
  int label_id;
  double cost;
  int nbOverlap;
} Triple;


bool decreaseCost( void *tl, void *tr )
{
  return ( reinterpret_cast< Triple * >( tl ) )->cost < ( reinterpret_cast< Triple * >( tr ) )->cost;
}

inline void actualizeTabuCandidateList( int m, int iteration, int nbOverlap, int *candidateListSize,
                                        double candidateBaseFactor, double *candidateFactor,
                                        int minCandidateListSize, double reductionFactor,
                                        int minTabuTSize, double tabuFactor, int *tenure, int n )
{

  if ( *candidateFactor > candidateBaseFactor )
    *candidateFactor /= reductionFactor;

  if ( iteration % m == 0 )
  {
    *tenure = minTabuTSize + static_cast< int >( tabuFactor * nbOverlap );
    *candidateListSize = minCandidateListSize + static_cast< int >( *candidateFactor * nbOverlap );

    if ( *candidateListSize > n )
      *candidateListSize = n;
  }

}


inline void actualizeCandidateList( int nbOverlap, int *candidateListSize, double,
                                    double *candidateFactor, int minCandidateListSize, double growingFactor, int n )
{
  if ( *candidateListSize < n )
    *candidateFactor = *candidateFactor * growingFactor;
  *candidateListSize = minCandidateListSize + static_cast< int >( *candidateFactor * nbOverlap );

  if ( *candidateListSize > n )
    *candidateListSize = n;
}




typedef struct
{
  LabelPosition *lp = nullptr;
  Triple **candidates = nullptr;
  double *labelPositionCost = nullptr;
  int *nbOlap = nullptr;
  double diff_cost;
  int *featWrap = nullptr;
  int *sol = nullptr;
  int borderSize;
} UpdateContext;

bool updateCandidatesCost( LabelPosition *lp, void *context )
{
  UpdateContext *ctx = reinterpret_cast< UpdateContext * >( context );

  if ( ctx->lp->isInConflict( lp ) )
  {
    ctx->labelPositionCost[lp->getId()] += ctx->diff_cost;
    if ( ctx->diff_cost > 0 )
      ctx->nbOlap[lp->getId()]++;
    else
      ctx->nbOlap[lp->getId()]--;

    int feat_id = ctx->featWrap[ctx->lp->getProblemFeatureId()];
    int feat_id2;
    if ( feat_id >= 0 && ctx->sol[feat_id] == lp->getId() ) // this label is in use
    {
      if ( ( feat_id2 = feat_id - ctx->borderSize ) >= 0 )
      {
        ctx->candidates[feat_id2]->cost += ctx->diff_cost;
        ctx->candidates[feat_id2]->nbOverlap--;
      }
    }
  }
  return true;
}




double Problem::popmusic_tabu( SubPart *part )
{
//+++++++++++++++++++gpl++++++++++++++++++++++++++++++++++
if(gplDebugger){
 cout<< "popmusic_tabu"<<endl;
}
//-------------------gpl-----------------------------------
  int probSize = part->probSize;
  int borderSize = part->borderSize;
  int subSize = part->subSize;
  int *sub = part->sub;
  int *sol = part->sol;

  Triple **candidateList = new Triple*[probSize];
  Triple **candidateListUnsorted = new Triple*[probSize];

  int *best_sol = new int[subSize];

  double cur_cost;
  double best_cost;
  double initial_cost;

  int *tabu_list = new int[probSize];

  int i;
  int j;

  int itwImp;
  int it = 0;
  int max_it;
  int stop_it;

  double delta;
  double delta_min;
  bool authorized;

  int feat_id;
  int feat_sub_id;
  int label_id;
  int p;

  int choosed_feat;
  int choosed_label;
  int candidateId;

  int nbOverlap = 0;
  //int nbOverlapLabel;


  int tenure = 10;  //
  int m = 50; // m   [10;50]

  int minTabuTSize = 9; // min_t [2;10]
  double tabuFactor = 0.5; // p_t [0.1;0.8]

  int minCandidateListSize = 18; // min_c   [2;20]

  double candidateBaseFactor = 0.73; // p_base  [0.1;0.8]

  double growingFactor = 15; // fa  [5;20]
  double reductionFactor = 1.3; // f_r [1.1;1.5]

  int candidateListSize = minCandidateListSize;
  double candidateFactor = candidateBaseFactor;


  int first_lp;

  //double EPSILON = 0.000001;
  max_it = probSize * pal->tabuMaxIt;
  itwImp = probSize * pal->tabuMinIt;
  stop_it = itwImp;

  cur_cost = 0.0;
  nbOverlap = 0;


  int lp;
  for ( i = 0; i < subSize; i++ )
    featWrap[sub[i]] = i;

  for ( i = 0; i < subSize; i++ )
  {
    j = featStartId[sub[i]];
    for ( lp = 0; lp < featNbLp[sub[i]]; lp++ )
    {
      it = j + lp;
      labelPositionCost[it] = compute_feature_cost( part, i, it, & ( nbOlap[it] ) );
    }
  }

  first_lp = ( displayAll ? 0 : -1 );
  for ( i = 0; i < probSize; i++ )
  {

    tabu_list[i] = -1; // nothing is tabu

    candidateList[i] = new Triple();
    candidateList[i]->feat_id = i + borderSize;
    candidateList[i]->label_id = sol[i + borderSize];

    candidateListUnsorted[i] = candidateList[i];

    if ( sol[i + borderSize] >= 0 )
    {
      j = sol[i + borderSize];
      candidateList[i]->cost = labelPositionCost[j];
      candidateList[i]->nbOverlap = nbOlap[j];
    }
    else
    {
      candidateList[i]->cost = inactiveCost[sub[i + borderSize]];
      candidateList[i]->nbOverlap = 1;
    }

  }


  for ( i = 0; i < subSize; i++ )
  {
    if ( sol[i] == -1 )
    {
      cur_cost += inactiveCost[sub[i]];
    }
    else
    {
      nbOverlap += nbOlap[sol[i]];
      cur_cost += labelPositionCost[sol[i]];
    }
  }

  Util::sort( reinterpret_cast< void ** >( candidateList ), probSize, decreaseCost );

  best_cost = cur_cost;
  initial_cost = cur_cost;
  memcpy( best_sol, sol, sizeof( int ) * ( subSize ) );

  // START TABU

  it = 0;
  while ( it < stop_it && best_cost >= EPSILON )
  {
    actualizeTabuCandidateList( m, it, nbOverlap, &candidateListSize, candidateBaseFactor, &candidateFactor, minCandidateListSize, reductionFactor, minTabuTSize, tabuFactor, &tenure, probSize );
    delta_min     = std::numeric_limits<double>::max();
    choosed_feat  = -1;
    choosed_label = -2;
    candidateId   = -1;

    // foreach retained Candidate
    for ( i = 0; i < candidateListSize; i++ )
    {
      feat_id     = candidateList[i]->feat_id;
      feat_sub_id = sub[feat_id];
      label_id    = candidateList[i]->label_id;

      int oldPos  = ( label_id < 0 ? -1 : label_id - featStartId[feat_sub_id] );


      p = featNbLp[feat_sub_id];

      /* label -1 means inactive feature */
      // foreach labelposition of feature minus the one in the solution
      for ( j = first_lp; j < p; j++ )
      {
        if ( j != oldPos )
        {

          if ( sol[feat_id] < 0 )
          {
            delta = -inactiveCost[feat_sub_id];
          }
          else
          {
            delta = -labelPositionCost[sol[feat_id]];
            delta -= nbOlap[sol[feat_id]] * ( inactiveCost[feat_sub_id] + mLabelPositions.at( label_id )->cost() );
          }

          if ( j >= 0 )
          {
            delta += labelPositionCost[featStartId[feat_sub_id] + j];
            delta += nbOlap[featStartId[feat_sub_id] + j] * ( inactiveCost[feat_sub_id] + mLabelPositions.at( featStartId[feat_sub_id] + j )->cost() );
          }
          else
          {
            delta += inactiveCost[feat_sub_id];
          }

          // move is authorized wether the feat isn't taboo or whether the move give a new best solution
          authorized = ( tabu_list[feat_id - borderSize] <= it ) || ( cur_cost + delta < best_cost );

          if ( delta < delta_min && authorized )
          {
            choosed_feat = feat_id;

            if ( j == -1 )
              choosed_label = -1;
            else
              choosed_label = featStartId[feat_sub_id] + j;

            delta_min = delta;
            candidateId = i;
          }
        }
      }
    }

    // if a modification has been retained
    if ( choosed_feat >= 0 )
    {
      // update the solution and update tabu list
      int old_label = sol[choosed_feat];

      tabu_list[choosed_feat - borderSize] = it + tenure;
      sol[choosed_feat] = choosed_label;
      candidateList[candidateId]->label_id = choosed_label;

      if ( old_label != -1 )
        mLabelPositions.at( old_label )->removeFromIndex( candidates_subsol );

      /* re-compute all labelpositioncost that overlap with old an new label */
      double local_inactive = inactiveCost[sub[choosed_feat]];

      if ( choosed_label != -1 )
      {
        candidateList[candidateId]->cost = labelPositionCost[choosed_label];
        candidateList[candidateId]->nbOverlap = nbOlap[choosed_label];
      }
      else
      {
        candidateList[candidateId]->cost = local_inactive;
        candidateList[candidateId]->nbOverlap = 1;
      }

      cur_cost += delta_min;

      double amin[2];
      double amax[2];
      LabelPosition *lp = nullptr;

      UpdateContext context;

      context.candidates = candidateListUnsorted;
      context.labelPositionCost = labelPositionCost;
      context.nbOlap = nbOlap;
      context.featWrap = featWrap;
      context.sol = sol;
      context.borderSize = borderSize;

      if ( old_label >= 0 )
      {
        lp = mLabelPositions.at( old_label );

        lp->getBoundingBox( amin, amax );

        context.diff_cost = -local_inactive - lp->cost();
        context.lp = lp;

        candidates->Search( amin, amax, updateCandidatesCost, &context );
      }

      if ( choosed_label >= 0 )
      {
        lp = mLabelPositions.at( choosed_label );

        lp->getBoundingBox( amin, amax );

        context.diff_cost = local_inactive + lp->cost();
        context.lp = lp;


        candidates->Search( amin, amax, updateCandidatesCost, &context );

        lp->insertIntoIndex( candidates_subsol );
      }

      Util::sort( reinterpret_cast< void ** >( candidateList ), probSize, decreaseCost );

      if ( best_cost - cur_cost > EPSILON ) // new best sol
      {
        best_cost = cur_cost;
        memcpy( best_sol, sol, sizeof( int ) * ( subSize ) );
        stop_it = it + itwImp;
        if ( stop_it > max_it )
          stop_it = max_it;
      }
    }
    else
    {
      /* no feature was selected : increase candidate list size*/
      actualizeCandidateList( nbOverlap, &candidateListSize, candidateBaseFactor,
                              &candidateFactor, minCandidateListSize, growingFactor, probSize );
    }
    it++;
  }

  memcpy( sol, best_sol, sizeof( int ) * ( subSize ) );

  for ( i = 0; i < subSize; i++ )
    featWrap[sub[i]] = -1;

  for ( i = 0; i < probSize; i++ )
    delete candidateList[i];

  delete[] candidateList;
  delete[] candidateListUnsorted;
  delete[] best_sol;
  delete[] tabu_list;

  /* Returns delta */
  return initial_cost - best_cost;
}





typedef struct
{
  LabelPosition *lp = nullptr;
  int *tmpsol = nullptr;
  int *featWrap = nullptr;
  int *feat = nullptr;
  int borderSize;
  QLinkedList<ElemTrans *> *currentChain;
  QLinkedList<int> *conflicts;
  double *delta_tmp = nullptr;
  double *inactiveCost = nullptr;

} ChainContext;


bool chainCallback( LabelPosition *lp, void *context )
{
  ChainContext *ctx = reinterpret_cast< ChainContext * >( context );

  if ( lp->isInConflict( ctx->lp ) )
  {
    int feat, rfeat;
    bool sub = nullptr != ctx->featWrap;

    feat = lp->getProblemFeatureId();
    if ( sub )
    {
      rfeat = feat;
      feat = ctx->featWrap[feat];
    }
    else
      rfeat = feat;

    if ( feat >= 0 && ctx->tmpsol[feat] == lp->getId() )
    {
      if ( sub && feat < ctx->borderSize )
      {
        throw - 2;
      }
    }

    // is there any cycles ?
    QLinkedList< ElemTrans * >::iterator cur;
    for ( cur = ctx->currentChain->begin(); cur != ctx->currentChain->end(); ++cur )
    {
      if ( ( *cur )->feat == feat )
      {
        throw - 1;
      }
    }

    if ( !ctx->conflicts->contains( feat ) )
    {
      ctx->conflicts->append( feat );
      *ctx->delta_tmp += lp->cost() + ctx->inactiveCost[rfeat];
    }
  }
  return true;
}

inline Chain *Problem::chain( SubPart *part, int seed )
{
//+++++++++++++++++++gpl++++++++++++++++++++++++++++
if(gplDebugger){
 std::cout<< "chain"<<endl;
}
//-------------------gpl----------------------------
  int i;
  int j;

  int lid;

  //int probSize   = part->probSize;
  int borderSize = part->borderSize;
  int subSize    = part->subSize;
  int *sub       = part->sub;
  int *sol       = part->sol;
  int subseed;

  double delta;
  double delta_min;
  double delta_best = std::numeric_limits<double>::max();
  double delta_tmp;

  int next_seed;
  int retainedLabel;
  Chain *retainedChain = nullptr;

  int max_degree = pal->ejChainDeg;

  int seedNbLp;

  QLinkedList<ElemTrans *> *currentChain = new QLinkedList<ElemTrans *>;
  QLinkedList<int> *conflicts = new QLinkedList<int>;

  int *tmpsol = new int[subSize];
  memcpy( tmpsol, sol, sizeof( int ) *subSize );

  LabelPosition *lp = nullptr;
  double amin[2];
  double amax[2];

  ChainContext context;
  context.featWrap = featWrap;
  context.borderSize = borderSize;
  context.tmpsol = tmpsol;
  context.inactiveCost = inactiveCost;
  context.feat = nullptr;
  context.currentChain = currentChain;
  context.conflicts = conflicts;
  context.delta_tmp = &delta_tmp;

  delta = 0;
  while ( seed != -1 )
  {
    subseed = sub[seed];
    seedNbLp = featNbLp[subseed];
    delta_min = std::numeric_limits<double>::max();
    next_seed = -1;
    retainedLabel = -2;


    if ( tmpsol[seed] == -1 )
      delta -= inactiveCost[subseed];
    else
      delta -= mLabelPositions.at( tmpsol[seed] )->cost();

    // TODO modify to handle displayAll param
    for ( i = -1; i < seedNbLp; i++ )
    {
      try
      {
        // Skip active label !
        if ( !( tmpsol[seed] == -1 && i == -1 ) && i + featStartId[subseed] != tmpsol[seed] )
        {
          if ( i != -1 )
          {
            lid = featStartId[subseed] + i;
            delta_tmp = delta;

            lp = mLabelPositions.at( lid );

            // evaluate conflicts graph in solution after moving seed's label
            lp->getBoundingBox( amin, amax );

            context.lp = lp;

            // search ative conflicts and count them
            candidates_subsol->Search( amin, amax, chainCallback, reinterpret_cast< void * >( &context ) );

            // no conflict -> end of chain
            if ( conflicts->isEmpty() )
            {
              if ( !retainedChain || delta + lp->cost() < delta_best )
              {

                if ( retainedChain )
                {
                  delete[] retainedChain->label;
                  delete[] retainedChain->feat;
                }
                else
                {
                  retainedChain = new Chain(); // HERE
                }

                delta_best = delta + lp->cost();

                retainedChain->degree = currentChain->size() + 1;
                retainedChain->feat  = new int[retainedChain->degree]; // HERE
                retainedChain->label = new int[retainedChain->degree]; // HERE
                QLinkedList< ElemTrans * >::iterator current = currentChain->begin();
                ElemTrans *move = nullptr;
                j = 0;
                while ( current != currentChain->end() )
                {
                  move = *current;
                  retainedChain->feat[j]  = move->feat;
                  retainedChain->label[j] = move->new_label;
                  ++current;
                  ++j;
                }
                retainedChain->feat[j] = seed;
                retainedChain->label[j] = lid;
                retainedChain->delta = delta + mLabelPositions.at( retainedChain->label[j] )->cost();
              }
            }

            // another feature can be ejected
            else if ( conflicts->size() == 1 )
            {
              if ( delta_tmp < delta_min )
              {
                delta_min = delta_tmp;
                retainedLabel = lid;
                next_seed = conflicts->takeFirst();
              }
              else
              {
                conflicts->takeFirst();
              }
            }
            else
            {

              // A lot of conflict : make them inactive and store chain
              Chain *newChain = new Chain();  // HERE
              newChain->degree = currentChain->size() + 1 + conflicts->size();
              newChain->feat  = new int[newChain->degree]; // HERE
              newChain->label = new int[newChain->degree]; // HERE
              QLinkedList<ElemTrans *>::iterator current = currentChain->begin();
              ElemTrans *move = nullptr;
              j = 0;
              while ( current != currentChain->end() )
              {
                move = *current;
                newChain->feat[j]  = move->feat;
                newChain->label[j] = move->new_label;
                ++current;
                ++j;
              }

              newChain->feat[j] = seed;
              newChain->label[j] = lid;
              newChain->delta = delta + mLabelPositions.at( newChain->label[j] )->cost();
              j++;


              while ( !conflicts->isEmpty() )
              {
                int ftid = conflicts->takeFirst();
                newChain->feat[j] = ftid;
                newChain->label[j] = -1;
                newChain->delta += inactiveCost[sub[ftid]];
                j++;
              }

              if ( newChain->delta < delta_best )
              {
                if ( retainedChain )
                  delete_chain( retainedChain );

                delta_best = newChain->delta;
                retainedChain = newChain;
              }
              else
              {
                delete_chain( newChain );
              }
            }
          }
          else   // Current label == -1   end of chain ...
          {
            if ( !retainedChain || delta + inactiveCost[subseed] < delta_best )
            {
              if ( retainedChain )
              {
                delete[] retainedChain->label;
                delete[] retainedChain->feat;
              }
              else
                retainedChain = new Chain(); // HERE

              delta_best = delta + inactiveCost[subseed];

              retainedChain->degree = currentChain->size() + 1;
              retainedChain->feat  = new int[retainedChain->degree]; // HERE
              retainedChain->label = new int[retainedChain->degree]; // HERE
              QLinkedList<ElemTrans *>::iterator current = currentChain->begin();
              ElemTrans *move = nullptr;
              j = 0;
              while ( current != currentChain->end() )
              {
                move = *current;
                retainedChain->feat[j]  = move->feat;
                retainedChain->label[j] = move->new_label;
                ++current;
                ++j;
              }
              retainedChain->feat[j] = seed;
              retainedChain->label[j] = -1;
              retainedChain->delta = delta + inactiveCost[subseed];
            }
          }
        }
      }
      catch ( int i )
      {
        Q_UNUSED( i );
        conflicts->clear();
      }
    } // end foreach labelposition

    if ( next_seed == -1 )
    {
      seed = -1;
    }
    else if ( currentChain->size() > max_degree )
    {
      seed = -1;
    }
    else
    {
      ElemTrans *et = new ElemTrans();
      et->feat  = seed;
      et->old_label = tmpsol[seed];
      et->new_label = retainedLabel;
      currentChain->append( et );

      if ( et->old_label != -1 )
      {
        mLabelPositions.at( et->old_label )->removeFromIndex( candidates_subsol );
      }

      if ( et->new_label != -1 )
      {
        mLabelPositions.at( et->new_label )->insertIntoIndex( candidates_subsol );
      }

      tmpsol[seed] = retainedLabel;
      delta += mLabelPositions.at( retainedLabel )->cost();
      seed = next_seed;
    }
  }

  while ( !currentChain->isEmpty() )
  {
    ElemTrans *et = currentChain->takeFirst();

    if ( et->new_label != -1 )
    {
      mLabelPositions.at( et->new_label )->removeFromIndex( candidates_subsol );
    }

    if ( et->old_label != -1 )
    {
      mLabelPositions.at( et->old_label )->insertIntoIndex( candidates_subsol );
    }

    delete et;
  }
  delete currentChain;

  delete[] tmpsol;
  delete conflicts;


  return retainedChain;
}


inline Chain *Problem::chain( int seed )
{
//+++++++++++++++++++gpl++++++++++++++++++++++++++++
if(gplDebugger){
 cout<<"chain"<<endl;
}
//-------------------gpl----------------------------
  int i;
  int j;

  int lid;

  double delta;
  double delta_min;
  double delta_best = std::numeric_limits<double>::max();
  double delta_tmp;

  int next_seed;
  int retainedLabel;
  Chain *retainedChain = nullptr;

  int max_degree = pal->ejChainDeg;

  int seedNbLp;

  QLinkedList<ElemTrans *> *currentChain = new QLinkedList<ElemTrans *>;
  QLinkedList<int> *conflicts = new QLinkedList<int>;

  int *tmpsol = new int[nbft];
  memcpy( tmpsol, sol->s, sizeof( int ) *nbft );

  LabelPosition *lp = nullptr;
  double amin[2];
  double amax[2];

  ChainContext context;
  context.featWrap = nullptr;
  context.borderSize = 0;
  context.tmpsol = tmpsol;
  context.inactiveCost = inactiveCost;
  context.feat = nullptr;
  context.currentChain = currentChain;
  context.conflicts = conflicts;
  context.delta_tmp = &delta_tmp;

  delta = 0;
  while ( seed != -1 )
  {
    seedNbLp = featNbLp[seed];
    delta_min = std::numeric_limits<double>::max();

    next_seed = -1;
    retainedLabel = -2;

    // sol[seed] is ejected
    if ( tmpsol[seed] == -1 )
      delta -= inactiveCost[seed];
    else
      delta -= mLabelPositions.at( tmpsol[seed] )->cost();

    for ( i = -1; i < seedNbLp; i++ )
    {
      try
      {
        // Skip active label !
        if ( !( tmpsol[seed] == -1 && i == -1 ) && i + featStartId[seed] != tmpsol[seed] )
        {
          if ( i != -1 ) // new_label
          {
            lid = featStartId[seed] + i;
            delta_tmp = delta;

            lp = mLabelPositions.at( lid );

            // evaluate conflicts graph in solution after moving seed's label
            lp->getBoundingBox( amin, amax );

            context.lp = lp;

            candidates_sol->Search( amin, amax, chainCallback, reinterpret_cast< void * >( &context ) );

            // no conflict -> end of chain
            if ( conflicts->isEmpty() )
            {
              if ( !retainedChain || delta + lp->cost() < delta_best )
              {
                if ( retainedChain )
                {
                  delete[] retainedChain->label;
                  delete[] retainedChain->feat;
                }
                else
                {
                  retainedChain = new Chain();
                }

                delta_best = delta + lp->cost();

                retainedChain->degree = currentChain->size() + 1;
                retainedChain->feat  = new int[retainedChain->degree];
                retainedChain->label = new int[retainedChain->degree];
                QLinkedList<ElemTrans *>::iterator current = currentChain->begin();
                ElemTrans *move = nullptr;
                j = 0;
                while ( current != currentChain->end() )
                {
                  move = *current;
                  retainedChain->feat[j]  = move->feat;
                  retainedChain->label[j] = move->new_label;
                  ++current;
                  ++j;
                }
                retainedChain->feat[j] = seed;
                retainedChain->label[j] = lid;
                retainedChain->delta = delta + lp->cost();
              }
            }

            // another feature can be ejected
            else if ( conflicts->size() == 1 )
            {
              if ( delta_tmp < delta_min )
              {
                delta_min = delta_tmp;
                retainedLabel = lid;
                next_seed = conflicts->takeFirst();
              }
              else
              {
                conflicts->takeFirst();
              }
            }
            else
            {

              // A lot of conflict : make them inactive and store chain
              Chain *newChain = new Chain();
              newChain->degree = currentChain->size() + 1 + conflicts->size();
              newChain->feat  = new int[newChain->degree];
              newChain->label = new int[newChain->degree];
              QLinkedList<ElemTrans *>::iterator current = currentChain->begin();
              ElemTrans *move = nullptr;
              j = 0;

              while ( current != currentChain->end() )
              {
                move = *current;
                newChain->feat[j]  = move->feat;
                newChain->label[j] = move->new_label;
                ++current;
                ++j;
              }

              // add the current candidates into the chain
              newChain->feat[j] = seed;
              newChain->label[j] = lid;
              newChain->delta = delta + mLabelPositions.at( newChain->label[j] )->cost();
              j++;

              // hide all conflictual candidates
              while ( !conflicts->isEmpty() )
              {
                int ftid = conflicts->takeFirst();
                newChain->feat[j] = ftid;
                newChain->label[j] = -1;
                newChain->delta += inactiveCost[ftid];
                j++;
              }

              if ( newChain->delta < delta_best )
              {
                if ( retainedChain )
                  delete_chain( retainedChain );

                delta_best = newChain->delta;
                retainedChain = newChain;
              }
              else
              {
                delete_chain( newChain );
              }
            }

          }
          else   // Current label == -1   end of chain ...
          {
            if ( !retainedChain || delta + inactiveCost[seed] < delta_best )
            {
              if ( retainedChain )
              {
                delete[] retainedChain->label;
                delete[] retainedChain->feat;
              }
              else
                retainedChain = new Chain();

              delta_best = delta + inactiveCost[seed];

              retainedChain->degree = currentChain->size() + 1;
              retainedChain->feat  = new int[retainedChain->degree];
              retainedChain->label = new int[retainedChain->degree];
              QLinkedList<ElemTrans *>::iterator current = currentChain->begin();
              ElemTrans *move = nullptr;
              j = 0;
              while ( current != currentChain->end() )
              {
                move = *current;
                retainedChain->feat[j]  = move->feat;
                retainedChain->label[j] = move->new_label;
                ++current;
                ++j;
              }
              retainedChain->feat[j] = seed;
              retainedChain->label[j] = -1;
              retainedChain->delta = delta + inactiveCost[seed];
            }
          }
        }
      }
      catch ( int i )
      {
        Q_UNUSED( i );
        conflicts->clear();
      }
    } // end foreach labelposition

    if ( next_seed == -1 )
    {
      seed = -1;
    }
    else if ( currentChain->size() > max_degree )
    {
      // Max degree reached
      seed = -1;
    }
    else
    {
      ElemTrans *et = new ElemTrans();
      et->feat  = seed;
      et->old_label = tmpsol[seed];
      et->new_label = retainedLabel;
      currentChain->append( et );

      if ( et->old_label != -1 )
      {
        mLabelPositions.at( et->old_label )->removeFromIndex( candidates_sol );
      }

      if ( et->new_label != -1 )
      {
        mLabelPositions.at( et->new_label )->insertIntoIndex( candidates_sol );
      }


      tmpsol[seed] = retainedLabel;
      delta += mLabelPositions.at( retainedLabel )->cost();
      seed = next_seed;
    }
  }


  while ( !currentChain->isEmpty() )
  {
    ElemTrans *et = currentChain->takeFirst();

    if ( et->new_label != -1 )
    {
      mLabelPositions.at( et->new_label )->removeFromIndex( candidates_sol );
    }

    if ( et->old_label != -1 )
    {
      mLabelPositions.at( et->old_label )->insertIntoIndex( candidates_sol );
    }

    delete et;
  }
  delete currentChain;

  delete[] tmpsol;
  delete conflicts;


  return retainedChain;
}

double Problem::popmusic_chain( SubPart *part )
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
if(gplDebugger){
 cout<< "popmusic_chain"<<endl;
}
//-------------------gpl-----------------------------------
  int i;
  //int j;

  int probSize   = part->probSize;
  int borderSize = part->borderSize;
  int subSize    = part->subSize;
  int *sub       = part->sub;
  int *sol       = part->sol;

  int *best_sol = new int[subSize];

  for ( i = 0; i < subSize; i++ )
  {
    featWrap[sub[i]] = i;
    best_sol[i] = sol[i];
  }

  double initial_cost;
  double cur_cost = 0;
  double best_cost = 0;

  // int nbOverlap = 0;

  int seed;

  int featOv;

  int lid;
  int fid;

  int *tabu_list = new int[subSize];

  Chain *current_chain = nullptr;

  //int itC;

  int it;
  int stop_it;
  int maxit;
  int itwimp; // iteration without improvment

  int tenure = pal->tenure;

  for ( i = 0; i < subSize; i++ )
  {
    cur_cost += compute_feature_cost( part, i, sol[i], &featOv );
    // nbOverlap += featOv;
  }

  initial_cost = cur_cost;
  best_cost = cur_cost;

  it = 0;

  maxit = probSize * pal->tabuMaxIt;

  itwimp = probSize * pal->tabuMinIt;

  stop_it = itwimp;

  // feature on border always are tabu
  for ( i = 0; i < borderSize; i++ )
    tabu_list[i] = maxit;   // border always are taboo

  for ( i = 0; i < probSize; i++ )
    tabu_list[i + borderSize] = -1; // others aren't

  while ( it < stop_it )
  {
    seed = ( it % probSize ) + borderSize;

    current_chain = chain( part, seed );
    if ( current_chain )
    {

      /* we accept a modification only if the seed is not tabu or
       * if the nmodification will generate a new best solution */
      if ( tabu_list[seed] < it || ( cur_cost + current_chain->delta ) - best_cost < 0.0 )
      {

        for ( i = 0; i < current_chain->degree; i++ )
        {
          fid = current_chain->feat[i];
          lid = current_chain->label[i];

          if ( sol[fid] >= 0 )
          {
            mLabelPositions.at( sol[fid] )->removeFromIndex( candidates_subsol );
          }
          sol[fid] = lid;

          if ( sol[fid] >= 0 )
          {
            mLabelPositions.at( lid )->insertIntoIndex( candidates_subsol );
          }

          tabu_list[fid] = it + tenure;
        }

        cur_cost += current_chain->delta;

        /* check if new solution is a new best solution */
        if ( best_cost - cur_cost > EPSILON )
        {
          best_cost = cur_cost;
          memcpy( best_sol, sol, sizeof( int ) *subSize );

          stop_it = ( it + itwimp > maxit ? maxit : it + itwimp );
        }
      }
      delete_chain( current_chain );
    }
    it++;
  }

  memcpy( sol, best_sol, sizeof( int ) *subSize );

  /*
  for (i=borderSize;i<subSize;i++){
     chain = chain (part, i);
     if (chain){
        if (chain->delta < 0.0){
           best_cost += chain->delta;
           for (j=0;j<chain->degree;j++){
              fid = chain->feat[j];
              lid = chain->label[j];
              sol[fid] = lid;
           }
        }

        delete_chain(chain);
     }
  }
  */

  for ( i = 0; i < subSize; i++ )
    featWrap[sub[i]] = -1;

  delete[] best_sol;
  delete[] tabu_list;


  return initial_cost - best_cost;
}

double Problem::popmusic_tabu_chain( SubPart *part )
{
//+++++++++++++++++++gpl++++++++++++++++++++++++++++++++++++++++
if(gplDebugger){
 cout<<"popmusic_tabu_chain"<<endl;
}
//-------------------gpl----------------------------------------
  int i;

  int probSize   = part->probSize;
  int borderSize = part->borderSize;
  int subSize    = part->subSize;
  int *sub       = part->sub;
  int *sol       = part->sol;

  int *best_sol = new int[subSize];

  for ( i = 0; i < subSize; i++ )
  {
    featWrap[sub[i]] = i;
  }

  double initial_cost;
  double cur_cost = 0;
  double best_cost = 0;

  // int nbOverlap = 0;

  int seed;

  int featOv;

  int lid;
  int fid;

  int *tabu_list = new int[subSize];

  Chain *retainedChain = nullptr;
  Chain *current_chain = nullptr;

  int itC;

  int it;
  int stop_it;
  int maxit;
  int itwimp;

  int tenure = pal->tenure;

  //int deltaIt = 0;

  Triple **candidates = new Triple*[probSize];
  Triple **candidatesUnsorted = new Triple*[probSize];

  for ( i = 0; i < subSize; i++ )
  {
    cur_cost += compute_feature_cost( part, i, sol[i], &featOv );
    // nbOverlap += featOv;
  }

  initial_cost = cur_cost;
  best_cost = cur_cost;

  it = 0;

  maxit = probSize * pal->tabuMaxIt;

  itwimp = probSize * pal->tabuMinIt;

  stop_it = itwimp;

  for ( i = 0; i < borderSize; i++ )
    tabu_list[i] = maxit;


  for ( i = 0; i < probSize; i++ )
  {
    tabu_list[i + borderSize] = -1;

    candidates[i] = new Triple();
    candidates[i]->feat_id = i + borderSize;
    candidatesUnsorted[i] = candidates[i];

    candidates[i]->cost = ( sol[i + borderSize] == -1 ? inactiveCost[i + borderSize] : mLabelPositions.at( sol[i + borderSize] )->cost() );
  }

  Util::sort( reinterpret_cast< void ** >( candidates ), probSize, decreaseCost );

  int candidateListSize;
  candidateListSize = int ( pal->candListSize * static_cast< double >( probSize ) + 0.5 );

  if ( candidateListSize > probSize )
    candidateListSize = probSize;

  while ( it < stop_it )
  {
    retainedChain = nullptr;

    for ( itC = 0; itC < candidateListSize; itC++ )
    {
      seed = candidates[itC]->feat_id;

      current_chain = chain( part, seed );

      if ( current_chain )
      {
        // seed is not taboo or chain give us a new best solution
        if ( tabu_list[seed] < it || ( cur_cost + current_chain->delta ) - best_cost < 0.0 )
        {
          if ( !retainedChain )
          {
            retainedChain = current_chain;
          }
          else if ( current_chain->delta - retainedChain->delta < EPSILON )
          {
            delete_chain( retainedChain );
            retainedChain = current_chain;
          }
          else
          {
            delete_chain( current_chain );
          }
        }
        else
        {
          delete_chain( current_chain );
        }
      }
    } // end foreach candidates

    if ( retainedChain /*&& retainedChain->delta <= 0*/ )
    {
      for ( i = 0; i < retainedChain->degree; i++ )
      {
        fid = retainedChain->feat[i];
        lid = retainedChain->label[i];

        if ( sol[fid] >= 0 )
          mLabelPositions.at( sol[fid] )->removeFromIndex( candidates_subsol );

        sol[fid] = lid;

        if ( lid >= 0 )
          mLabelPositions.at( lid )->insertIntoIndex( candidates_subsol );

        tabu_list[fid] = it + tenure;
        candidatesUnsorted[fid - borderSize]->cost = ( lid == -1 ? inactiveCost[sub[fid]] : mLabelPositions.at( lid )->cost() );

      }

      cur_cost += retainedChain->delta;

      delete_chain( retainedChain );

      if ( best_cost - cur_cost > EPSILON )
      {
        best_cost = cur_cost;
        memcpy( best_sol, sol, sizeof( int ) * subSize );

        stop_it = ( it + itwimp > maxit ? maxit : it + itwimp );
      }
      Util::sort( reinterpret_cast< void ** >( candidates ), probSize, decreaseCost );
    }
    it++;
  }

  memcpy( sol, best_sol, sizeof( int ) *subSize );

  for ( i = 0; i < probSize; i++ )
    delete candidates[i];

  delete[] candidates;
  delete[] candidatesUnsorted;

  for ( i = 0; i < subSize; i++ )
    featWrap[sub[i]] = -1;

  delete[] best_sol;
  delete[] tabu_list;


  return initial_cost - best_cost;
}

bool checkCallback( LabelPosition *lp, void *ctx )
{
  QLinkedList<LabelPosition *> *list = reinterpret_cast< QLinkedList<LabelPosition *>* >( ctx );
  list->append( lp );

  return true;
}


void Problem::check_solution()
{
  int *solution = new int[nbft];

  double amin[2];
  double amax[2];

  amin[0] = bbox[0];
  amin[1] = bbox[1];
  amax[0] = bbox[2];
  amax[1] = bbox[3];

  QLinkedList<LabelPosition *> *list = new QLinkedList<LabelPosition *>;

  candidates_sol->Search( amin, amax, checkCallback, reinterpret_cast< void * >( list ) );

  int i;
  int nbActive = 0;
  for ( i = 0; i < nbft; i++ )
  {
    solution[i] = -1;
    if ( sol->s[i] >= 0 )
      nbActive++;
  }

  if ( list->size() != nbActive )
  {
    // Error in solution !!!!
  }

  while ( !list->isEmpty() )
  {
    LabelPosition *lp = list->takeFirst();
    int probFeatId = lp->getProblemFeatureId();

    solution[probFeatId] = lp->getId();
  }

  delete [] solution;
}

typedef struct _nokContext
{
  LabelPosition *lp = nullptr;
  bool *ok = nullptr;
  int *wrap = nullptr;
} NokContext;

bool nokCallback( LabelPosition *lp, void *context )
{
  NokContext *ctx = reinterpret_cast< NokContext *>( context );
  LabelPosition *lp2 = ctx->lp;
  bool *ok = ctx->ok;
  int *wrap = ctx->wrap;

  if ( lp2->isInConflict( lp ) )
  {
    if ( wrap )
    {
      ok[wrap[lp->getProblemFeatureId()]] = false;
    }
    else
    {
      ok[lp->getProblemFeatureId()] = false;
    }
  }

  return true;
}

void Problem::chain_search(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev)
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++
auto start = std::chrono::system_clock::now();
if(gplDebugger){
 std::cout<<"chain_search"<<endl;
}
//-------------------gpl----------------------------------
  if ( nbft == 0 )
    return;

  int i;
  int seed;
  bool *ok = new bool[nbft];
  int fid;
  int lid;
  int popit = 0;
  double amin[2];
  double amax[2];

  NokContext context;
  context.ok = ok;
  context.wrap = nullptr;

  Chain *retainedChain = nullptr;

  featWrap = nullptr;

  std::fill( ok, ok + nbft, false );

  //initialization();
  init_sol_falp(performance,false,is_initial,my_solution_prev);

  //check_solution();
  solution_cost();

  int iter = 0;

  while ( true )
  {

    //check_solution();

    for ( seed = ( iter + 1 ) % nbft;
          ok[seed] && seed != iter;
          seed = ( seed + 1 ) % nbft )
      ;

    // All seeds are OK
    if ( seed == iter )
    {
      break;
    }

    iter = ( iter + 1 ) % nbft;
    retainedChain = chain( seed );

    if ( retainedChain && retainedChain->delta < - EPSILON )
    {
      // apply modification
      for ( i = 0; i < retainedChain->degree; i++ )
      {
        fid = retainedChain->feat[i];
        lid = retainedChain->label[i];

        if ( sol->s[fid] >= 0 )
        {
          LabelPosition *old = mLabelPositions.at( sol->s[fid] );
          old->removeFromIndex( candidates_sol );

          old->getBoundingBox( amin, amax );

          context.lp = old;
          candidates->Search( amin, amax, nokCallback, &context );
        }

        sol->s[fid] = lid;

        if ( sol->s[fid] >= 0 )
        {
          mLabelPositions.at( lid )->insertIntoIndex( candidates_sol );
        }

        ok[fid] = false;
      }
      sol->cost += retainedChain->delta;
    }
    else
    {
      // no chain or the one is not god enough
      ok[seed] = true;
    }

    delete_chain( retainedChain );
    popit++;
  }

  solution_cost();
  performance.solutionWeight = sol->cost;
  delete[] ok;
  int solutionCount = 0;
  for(int i = 0; i < nbft; i++){
    if(sol->s[i]  != -1) solutionCount++;
  }
  performance.name = "chain_search";
  performance.solutionSize = solutionCount;
  auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
  bool testini = is_initial;
performance.remainingLabels = cacheSolution(true,is_initial,my_solution_prev);
if(testini) assert(performance.remainingLabels == 0);
  if(testPrinter){
    performance.print();
  }
}

bool Problem::compareLabelArea( pal::LabelPosition *l1, pal::LabelPosition *l2 )
{
  return l1->getWidth() * l1->getHeight() > l2->getWidth() * l2->getHeight();
}

QList<LabelPosition *> Problem::getSolution( bool returnInactive )
{
  int i;
  QList<LabelPosition *> solList;

  if ( nbft == 0 )
  {
    return solList;
  }

  for ( i = 0; i < nbft; i++ )
  {
    if ( sol->s[i] != -1 )
    {
      solList.push_back( mLabelPositions.at( sol->s[i] ) ); // active labels
    }
    else if ( returnInactive
              || mLabelPositions.at( featStartId[i] )->getFeaturePart()->layer()->displayAll()
              || mLabelPositions.at( featStartId[i] )->getFeaturePart()->alwaysShow() )
    {
      solList.push_back( mLabelPositions.at( featStartId[i] ) ); // unplaced label
    }
  }

  // if features collide, order by size, so smaller ones appear on top
  if ( returnInactive )
  {
    std::sort( solList.begin(), solList.end(), compareLabelArea );
  }

  return solList;
}

PalStat *Problem::getStats()
{
  int i, j;

  PalStat *stats = new PalStat();

  stats->nbObjects = nbft;
  stats->nbLabelledObjects = 0;

  stats->nbLayers = nbLabelledLayers;
  stats->layersNbObjects = new int[stats->nbLayers];
  stats->layersNbLabelledObjects = new int[stats->nbLayers];

  for ( i = 0; i < stats->nbLayers; i++ )
  {
    stats->layersName << labelledLayersName[i];
    stats->layersNbObjects[i] = 0;
    stats->layersNbLabelledObjects[i] = 0;
  }

  QString lyrName;
  int k;
  for ( i = 0; i < nbft; i++ )
  {
    lyrName = mLabelPositions.at( featStartId[i] )->getFeaturePart()->feature()->provider()->name();
    k = -1;
    for ( j = 0; j < stats->nbLayers; j++ )
    {
      if ( lyrName == stats->layersName[j] )
      {
        k = j;
        break;
      }
    }
    if ( k != -1 )
    {
      stats->layersNbObjects[k]++;
      if ( sol->s[i] >= 0 )
      {
        stats->layersNbLabelledObjects[k]++;
        stats->nbLabelledObjects++;
      }
    }
    else
    {
      // unknown layer
    }
  }

  return stats;
}

void Problem::solution_cost()
{
  sol->cost = 0.0;
  nbActive = 0;

  int nbOv;

  int i;

  LabelPosition::CountContext context;
  context.inactiveCost = inactiveCost;
  context.nbOv = &nbOv;
  context.cost = &sol->cost;
  double amin[2];
  double amax[2];
  LabelPosition *lp = nullptr;

  int nbHidden = 0;

  for ( i = 0; i < nbft; i++ )
  {
    if ( sol->s[i] == -1 )
    {
      sol->cost += inactiveCost[i];
      nbHidden++;
    }
    else
    {
      nbOv = 0;
      lp = mLabelPositions.at( sol->s[i] );

      lp->getBoundingBox( amin, amax );

      context.lp = lp;
      candidates_sol->Search( amin, amax, LabelPosition::countFullOverlapCallback, &context );

      sol->cost += lp->cost();

      if ( nbOv == 0 )
        nbActive++;
    }
  }
}
//+++++++++++++++++++++++++++gpl-algorithms++++simple++++++++++++++++++++++++++++++++++++++++++++++
// another implementation of simple
// no test yet
 /*typedef struct
 {
   bool* labelList = nullptr;
   LabelPosition *lp = nullptr;
   RTree <LabelPosition *, double, 2, double> *candidates;
 } SimpleContext;

bool simpleConflictCallback(LabelPosition *lp, void *ctx){
 SimpleContext *context = reinterpret_cast< SimpleContext * >( ctx );
 LabelPosition *choosenLp = context->lp;
 bool* list = context->labelList; 
 int id = lp->getId(); 
 if(id > choosenLp->getId() && choosenLp->isInConflict(lp)){
    list[id] = false;
 }
  return true;
}
void Problem::simple()
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
    if(gplDebugger){
        QgsLogger::QgsDebugMsg( QStringLiteral( "simple"));
    }
//-------------------gpl-----------------------------------
    int i, j;
    int label;
    double amin[2];
    double amax[2];
    QVector<bool> labelList(all_nblp,true);
    LabelPosition *lp = nullptr;
    SimpleContext *context = new SimpleContext();
    context->candidates = candidates;
    context->labelList = labelList.data();
    bool setF = false;
    init_sol_empty();
     for ( i = 0; i < nbft; i++){
        setF= false;
        for ( j = 0; j < featNbLp[i]; j++ )
        {
          if(setF) break;Simple
          label = featStaSimple
            if(labelList[Simple
            lp=mLabelPosiSimple
            if ( lp->getId() != label )
            {
                std::cerr << "simple wrong";
            }
            int probFeatId = lp->getProblemFeatureId();
            sol->s[probFeatId] = label; 
            lp->getBoundingBox( amin, amax); 
            //ignore all its overlapps;
            context->lp = lp;
            setF= true;
            candidates->Search( amin, amax, simpleConflictCallback ,reinterpret_cast< void * >( context ));        
        }simple
    }
   if ( dsimple
   {
     int nbOverlap;
     int start_p;
     LabelPosition *retainedLabel = nullptr;
     int p;
 
     for ( i = 0; i < nbft; i++ ) // forearch hidden feature
     {
       if ( sol->s[i] == -1 )
       {
         nbOverlap = std::numeric_limits<int>::max();
         start_p = featStartId[i];
         for ( p = 0; p < featNbLp[i]; p++ )
         {
           lp = mLabelPositions.at( start_p + p );
           lp->resetNumOverlaps();
 
           lp->getBoundingBox( amin, amax );
 
 
           candidates_sol->Search( amin, amax, LabelPosition::countOverlapCallback, lp );
 
           if ( lp->getNumOverlaps() < nbOverlap )
           {
             retainedLabel = lp;
             nbOverlap = lp->getNumOverlaps();
           }
         }
         sol->s[i] = retainedLabel->getId();
 
         retainedLabel->insertIntoIndex( candidates_sol );
 
       }
     }
   }
}*/
//------------------------gpl-algorithms---------------simple--------------------------------------
//+++++++++++++++++++++++++++gpl-algorithms+++simple++++++++++++++++++++++++++++++++++++++++++++
typedef struct
 {
   QSet<int> labelList;
   LabelPosition *lp = nullptr;
   RTree <LabelPosition *, double, 2, double> *candidates;
 } SimpleContext;

bool simpleConflictCallback(LabelPosition *lp, void *ctx){
 SimpleContext *context = reinterpret_cast< SimpleContext * >( ctx );
 LabelPosition *choosenLp = context->lp; 
 int id = lp->getId(); 
 if(id > choosenLp->getId()&&choosenLp->isInConflict(lp)){
    context->labelList.insert(id);
 }
  return true;
}
// greedy algorithm 
// no graph representation used
// the order of feature counts
void Problem::simple(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev)
{
//+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
    if(gplDebugger){
       cout<< "simple"<<endl;
    }
    int solutionCount= 0;
    if(gplPrinter){
      printCost();
    }
    auto start = std::chrono::system_clock::now();
//-------------------gpl-----------------------------------
    int i, j;
    int label;
    double amin[2];
    double amax[2];
    bool setF = false;
    LabelPosition *lp = nullptr;
    SimpleContext *context = new SimpleContext();
    context->candidates = candidates;
    init_sol_empty();
     for ( i = 0; i < nbft; i++){
        int offset = getCached(i,my_solution_prev);
        if(offset >-1){
          label = featStartId[i]+offset;
          if(!context->labelList.contains(label)){
            lp=mLabelPositions.at(label);
            if ( lp->getId() != label )
            {
                std::cerr << "simple wrong";
            }
            int probFeatId = lp->getProblemFeatureId();
            sol->s[probFeatId] = label; 
            lp->getBoundingBox( amin, amax); 
            //ignore all its overlapps;
            context->lp = lp;
            setF = true;
            candidates->Search( amin, amax, simpleConflictCallback ,reinterpret_cast< void * >( context ));   
            continue;
          }
        }
        setF= false;
        // initial solution 
        for ( j = 0; j < featNbLp[i]; j++ )
        {
            if(setF){
                break;
            }
            label = featStartId[i] + j;  
            if(context->labelList.contains(label)){
                 continue;
            }
            lp=mLabelPositions.at(label);
            if ( lp->getId() != label )
            {
                std::cerr << "simple wrong";
            }
            int probFeatId = lp->getProblemFeatureId();
            sol->s[probFeatId] = label; 
            solutionCount++;
            lp->getBoundingBox( amin, amax); 
            //ignore all its overlapps;
            context->lp = lp;
            setF = true;
            candidates->Search( amin, amax, simpleConflictCallback ,reinterpret_cast< void * >( context ));     
        }
    }
   //checkQgsfeatureID();
   delete context;
//+++++++++++++++++++gpl++++++++++++++++++++++++++++++++++
if(gplDebugger){
  solution_cost();
  cout<< sol->cost<<endl;
}
solution_cost();
performance.name = "simple";
performance.solutionWeight = sol->cost;
performance.solutionSize = solutionCount;
 bool testini = is_initial;
performance.remainingLabels = cacheSolution(true,is_initial,my_solution_prev);
if(testini) assert(performance.remainingLabels == 0);
auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
if(testPrinter){
  performance.print();
}
//-------------------gpl-----------------------------------
   if ( displayAll )
   {
     int nbOverlap;
     int start_p;
     LabelPosition *retainedLabel = nullptr;
     int p;
 
     for ( i = 0; i < nbft; i++ ) // forearch hidden feature
     {
       if ( sol->s[i] == -1 )
       {
         nbOverlap = std::numeric_limits<int>::max();
         start_p = featStartId[i];
         for ( p = 0; p < featNbLp[i]; p++ )
         {
           lp = mLabelPositions.at( start_p + p );
           lp->resetNumOverlaps();
 
           lp->getBoundingBox( amin, amax );
           candidates_sol->Search( amin, amax, LabelPosition::countOverlapCallback, lp );
 
           if ( lp->getNumOverlaps() < nbOverlap )
           {
             retainedLabel = lp;
             nbOverlap = lp->getNumOverlaps();
           }
         }
         sol->s[i] = retainedLabel->getId();
 
         retainedLabel->insertIntoIndex( candidates_sol );
 
       }
     }
   }
}
//------------------------gpl-algorithms-----------------simple---------------------------

//++++++++++++++++++++++++gpl-dataSruct++++set conflict graph+++++++++++++++++++++++++++++
typedef struct{
  Graph* conflictGraph;
  LabelPosition* lp;
  int* lpID= nullptr;
}CONFLICTContext;
// add bidirected edge for one conflict in conflictgraph
bool conflictCallBack( LabelPosition *lp, void *ctx ){
  CONFLICTContext* context = reinterpret_cast< CONFLICTContext* >( ctx );
  LabelPosition *lp2 = context->lp;
  int id = lp->getId();
  int id2 = *(context->lpID);
  if ( id > id2)
  {
    if(lp2->getProblemFeatureId()  == lp->getProblemFeatureId() || lp2->isInConflict( lp)){ 
      context->conflictGraph->addEdge(id,id2);
    }
  }
  return true;
}
// set conflictgraph according to the Rtree
//TODO: get the graph in the process of setting Rtree (in Rtree, it has counted once the number of overlappings)
void Problem::setConflictGraph(const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  conflictGraph = new Graph(nblp,all_nblp);
  //conflictGraph = new Graph(all_nblp);
  int i, j;
  int label;
  double amin[2];
  double amax[2];
  int lpID;
  int startID, endID;
  double weight;
  LabelPosition *lp = nullptr;
  CONFLICTContext * context = new CONFLICTContext();
  context-> lpID = &lpID;
  context->conflictGraph = conflictGraph;
  for ( i = 0; i < nbft; i++ ){
    startID = featStartId[i];
    for(j = 0; j < featNbLp[i]; j++){
        weight = inactiveCost[i]-mLabelPositions.at(startID+j)->cost();
        conflictGraph->addVertex(startID+j, weight);
    }
  }
  for ( i = 0; i < nbft; i++ ){
    startID = featStartId[i];
    endID = startID+ featNbLp[i];
    for ( j = startID; j <endID; j++ ){
      label = j;
      //add neighbors with smaller indices 
      for(int k = startID; k < j; k++){
        conflictGraph->addEdge(k,j);
        conflictGraph->addEdge(j,k);
      }
       lp = mLabelPositions.at( label );
       lpID = lp->getId();
       if ( lpID != label )
        {
          std::cerr << "mis wrong";
        }
       lp->getBoundingBox(amin, amax);
       context-> lp = lp;
       candidates->Search( amin, amax, conflictCallBack, reinterpret_cast< void * >(context) );   
    }
  }
  if(!is_initial){
    // build previou solution
    int qgsID;
    int offset;
    for ( i = 0; i < nbft; i++ ){
      startID = featStartId[i];
      qgsID = mLabelPositions.at(startID)->getFeaturePart()->feature()->id();
      std::unordered_map<int,int>::const_iterator got = my_solution_prev.find(qgsID);
      if (got != my_solution_prev.end()){
        offset = got->second;
        if(offset< featNbLp[i]){
          conflictGraph->addCache(featStartId[i] + offset);
        }
      }
    }
    conflictGraph->adjustWeights();
  }
}
// debug callback 
bool graphDebugCallBack( LabelPosition *lp, void *ctx ){
  CONFLICTContext* context = reinterpret_cast< CONFLICTContext* >( ctx );
  LabelPosition *lp2 = context->lp;
  int id = lp->getId();
  int id2 = *(context->lpID);
  if ( id != id2)
  {
    if(lp2->isInConflict( lp)){ 
      context->conflictGraph->containEdge(id,id2);
      context->conflictGraph->containEdge(id2,id);
    }
  }
  return true;

}
// debug code 
// check if all candidates of one feature are in conflict
//check if all and only the overlapped labels in RTree are connected in conflictgraph.
void Problem::debugConflictGraph(){
  conflictGraph->debugGraph();
  if(gplPrinter){
    conflictGraph->printGraph();
  }
  int i,j;
  int lpID;
  int start,end;
  LabelPosition *lp = nullptr;
  LabelPosition *lp2 = nullptr;
  int label,label2;
  vector<int> neighbors;
  for ( i = 0; i < nbft; i++ ){
    start = featStartId[i];
    end = start+featNbLp[i];
    for (j = start; j < end; j++)
    {
      label = j;
      assert(conflictGraph->containVertex(label));
      lp = mLabelPositions.at( label );
      lpID = lp->getId();
      if ( lpID != label )
      {
        std::cerr << "conflictGraph wrong";
      }
      for ( int p = 0; p < nbft; p++ ){
        for (int q = 0; q < featNbLp[p]; q++){
          label2 = featStartId[p]+q;
          lp2 = mLabelPositions.at( label2 );
          bool contains = conflictGraph->containEdge_label(label, label2);
          if(label == label2){
            assert(!contains);
            continue;
          } 
          bool conflict = (lp->isInConflict( lp2));
          bool neighbor = (lp->getFeaturePart() == lp2->getFeaturePart());
          assert((conflict || neighbor) == contains);
        }
      }
    }
  }
}
// bebug code for all MIS-getting-algorithms (set conflict graph and get maximal independet set)
//check the independency of labels 
// TODO: check indepedency directly in the original problem
void Problem::debugIndepdency( vector<int>& MIS){
  for(const auto &p: MIS){
    for(const auto &q: MIS){
      if(p == q) continue;
      assert(!conflictGraph->containEdge_label(p,q));
      if(mLabelPositions[p]->getFeaturePart() == mLabelPositions[q]->getFeaturePart()){
        std::cout<< "p: "<< p << endl;
         std::cout<< "q: "<< q<< endl;
        conflictGraph->printGraph();
      }
      assert(mLabelPositions[p]->getFeaturePart() != mLabelPositions[q]->getFeaturePart());
      assert(!mLabelPositions[p]->isInConflict( mLabelPositions[q] ));

    }
  }
}
//-------------------gpl-dataSruct----- set conflict graph -------------------------------------

//+++++++++++++++++++++++++++gpl-algorithms++++++++++++MIS++++++++++++++++++++++++++++++++++++++

void Problem::mis(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  //+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
    if(gplDebugger){
        cout<< "mis"<<endl;
    }
    auto start = std::chrono::system_clock::now();
//-------------------gpl-----------------------------------
  int label;
  int i,j;
  init_sol_empty();
  setConflictGraph(is_initial,my_solution_prev);
  if(gplDebugger){
    debugConflictGraph();
  }
  unordered_set<int> vertexCover = conflictGraph->getVertexCover_weighted(nblp, all_nblp);
  //unordered_set<int> vertexCover = conflictGraph->getVertexCover(nblp, all_nblp);
  vector<int> MIS;
  for ( i = 0; i < nbft; i++ ){
    for (j = 0; j < featNbLp[i]; j++ ){
      label = featStartId[i] + j;
      if(vertexCover.find(label) != vertexCover.end()) continue;
      MIS.push_back(label);
    }
  }
  if(gplDebugger){
    debugIndepdency(MIS);
  }
  setSolution(MIS);
  solution_cost();
  performance.solutionWeight = sol->cost;
  int remainingCount = cacheSolution(true,is_initial,my_solution_prev);
  int solutionCount = 0;
  for(int i = 0; i < nbft; i++){
    if(sol->s[i] != -1) solutionCount++;
  }
  performance.name = "mis";
  performance.solutionSize = solutionCount;
   bool testini = is_initial;
  performance.remainingLabels = remainingCount;
  if(testini) assert(performance.remainingLabels == 0);
  auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
  if(testPrinter){
    performance.print();
  }
}
// set the solution found using MIS-like algorithm
// a set of independent labels 
void Problem::setSolution(vector<int>& MIS){
  int i;
  double amin[2];
  double amax[2];
  LabelPosition *lp = nullptr;
  for(const auto &label: MIS){
    lp= mLabelPositions.at(label);
    if ( lp->getId() != label )
    {
      std::cerr << "setSolution wrong";
    }
    int probFeatId = lp->getProblemFeatureId();
    sol->s[probFeatId] = label; 
  }
  if ( displayAll )
  {
    int nbOverlap;
    int start_p;
    LabelPosition *retainedLabel = nullptr;
    int p;

    for ( i = 0; i < nbft; i++ ) // forearch hidden feature
    {
      if ( sol->s[i] == -1 )
      {
        nbOverlap = std::numeric_limits<int>::max();
        start_p = featStartId[i];
        for ( p = 0; p < featNbLp[i]; p++ )
        {
          lp = mLabelPositions.at( start_p + p );
          lp->resetNumOverlaps();

          lp->getBoundingBox( amin, amax );


          candidates_sol->Search( amin, amax, LabelPosition::countOverlapCallback, lp );

          if ( lp->getNumOverlaps() < nbOverlap )
          {
            retainedLabel = lp;
            nbOverlap = lp->getNumOverlaps();
          }
        }
        sol->s[i] = retainedLabel->getId();

        retainedLabel->insertIntoIndex( candidates_sol );

      }
    }
  }
}
// cover all edges 

//---------------------gpl-algorithms-----------------mis--------------------------


//+++++++++++++++++++++++++++gpl-algorithms++++++++++++MaxHS++++++++++++++++++++++++++++++++++++++
//get a maximum independent set using maxHS
void Problem::maxHS(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  //+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
    if(gplDebugger){
        cout<<"maxHS"<<endl;
    }
    auto start = std::chrono::system_clock::now();
//-------------------gpl-----------------------------------
  init_sol_empty();
  setConflictGraph(is_initial,my_solution_prev);
  if(gplDebugger){
    debugConflictGraph();
  }
  vector<int> KAMIS;
  conflictGraph->getMAXHS(KAMIS);
  if(gplDebugger){
    debugIndepdency(KAMIS);
  }
  setSolution(KAMIS);
  int remainingCount = cacheSolution(true,is_initial,my_solution_prev);
  int solutionCount = 0;
  for(int i = 0; i < nbft; i++){
    if(sol->s[i] != -1) solutionCount++;
  }
  solution_cost();
  performance.name = "maxHS";
  performance.solutionWeight = sol->cost;
  performance.solutionSize = solutionCount;
   bool testini = is_initial;
  performance.remainingLabels = remainingCount;
  if(testini) assert(performance.remainingLabels == 0);
  auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
    if(testPrinter){
    performance.print();
  }
}
//---------------------gpl-algorithms-----------------MaxHs---------------------------------------

//+++++++++++++++++++++++++++gpl-algorithms++++++++++++KAMIS++++++++++++++++++++++++++++++++++++++
//using KAMIS to solve it (maximum Indepdend set)
//currently as extern libery( file to file communication) 
// Unweighted
// Unstable
void Problem::kamis(test::Performance& performance,const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  //+++++++++++++++++++gpl+++++++++++++++++++++++++++++++++++
    if(gplDebugger){
        cout<< "kamis"<<endl;
    }
    auto start = std::chrono::system_clock::now();
//-------------------gpl-----------------------------------
  init_sol_empty();
  setConflictGraph(is_initial,my_solution_prev);
  if(gplDebugger){
    debugConflictGraph();
  }
  vector<int> KAMIS;
  conflictGraph->getKAMIS(KAMIS);
  if(gplDebugger){
    debugIndepdency(KAMIS);
  }
  setSolution(KAMIS);

  int offset;
  QgsFeatureId id;
  std::unordered_map<int,int>::const_iterator got;
  // count remaing labels
  int remainingCount = 0;
  if(!is_initial){
    for (int  i = 0; i < nbft; i++ ){
      offset = sol->s[i]-featStartId[i];
      if(offset >= 0){
        id = mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id();
        got = my_solution_prev.find (id);
        if ( got != my_solution_prev.end()){
          if(offset == got->second) remainingCount++;
        }
      }
    }
  }
  int solutionCount = 0;
  for(int i = 0; i < nbft; i++){
    if(sol->s[i] != -1) solutionCount++;
  }
  solution_cost();
  performance.name = "kamis";
  performance.solutionWeight = sol->cost;
  performance.solutionSize = solutionCount;
   bool testini = is_initial;
  performance.remainingLabels = remainingCount;
  if(testini) assert(performance.remainingLabels == 0);
  auto end = std::chrono::system_clock::now();
const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
 performance.time = duration;
  if(testPrinter){
    performance.print();
  }
}
//---------------------gpl-algorithms-----------------KAMIS---------------------------------------

//+++++++++++++++++++++++++++++modification+++++++++++++++++++++++++++++++++++++++++++++++++++++++
// cache the solution in solution_prev (map of qgsFeatureID and offset of the labelID)
int Problem:: cacheSolution(bool final,const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  if(!final) return -1;
  int offset;
  QgsFeatureId id;
  std::unordered_map<int,int>::const_iterator got;
  // count remaing labels
  int remainingCount = 0;
  if(!is_initial){
    for (int  i = 0; i < nbft; i++ ){
      offset = sol->s[i]-featStartId[i];
      if(offset >= 0){
        id = mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id();
        got = my_solution_prev.find (id);
        if ( got != my_solution_prev.end()){
          if(offset == got->second) remainingCount++;
        }
      }
    }
  }
  // store new solution
  my_solution_prev.clear();
  for (int  i = 0; i < nbft; i++ ){
    offset = sol->s[i]-featStartId[i];
    if(offset >= 0){
      // there are even features with same uniqueID
      /*std::unordered_map<int,int>::const_iterator got = solution_prev.find(mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id());
      assert(got == solution_prev.end());
      */
      my_solution_prev.insert({mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id(),offset});
    }
  }
  if(gplPrinter){
    cout<< my_solution_prev.size() << " are cached from the previous solution"<< endl;
    for (const auto &p : my_solution_prev) {
      std::cout << p.first << " => "<< p.second << endl;
    }
  }
  return remainingCount;
}
// get the offset to the startID (>= 0)
// if no cached avalid, return -1
int Problem:: getCached(int feat,unordered_map<int, int>& my_solution_prev){
  int offset;
  int qgsID = mLabelPositions.at(featStartId[feat])->getFeaturePart()->feature()->id();
  std::unordered_map<int,int>::const_iterator got = my_solution_prev.find(qgsID);
  if (got == my_solution_prev.end())
    return -1;
  offset = my_solution_prev[qgsID];
  if(offset< featNbLp[feat]){
    return offset;
  } 
  return -1;
}
//-----------------------------modification-------------------------------------------------------

//++++++++++++++++++++++++++++++++++Assertions for understanding QGIS+++++++++++++++++++++++++++++++++++
//pal::labelposition(for labeling), pal::featurePart(reduldent information TODO: integrate PAL with QGIS feature directly) and also Qgsfeature(assume has unique id).
//assume the id of qgsfeature is unique and static in same bounding map (different problem instance)
//QSet<int> QgsFeatureIDS_prev;
//int numF_prev;
//int numL_prev;
//bool init;
//unordered_map<int, int> solution_prev; 
//unordered_map<int, unordered_set<int>> candidates_prev;
void print(QSet<int>& prev){
  QSet<int>::const_iterator i;
    for (i = prev.begin(); i != prev.end(); ++i){
      cout<< *i<< " ";
    }
    cout<< endl;

}
void printLabelPosition(LabelPosition* lp){
  cout<< "label with id " << lp->getId()<< endl;
  cout<< "belongs to feature "<< lp->getFeaturePart()->feature()->id()<< endl;
 // cout<< "its coordinats are "<< lp->getX() << ", "<< lp->getY()<< endl;
  //cout<< "width and height are "<< lp->getWidth() << "." << lp->getHeight() << endl;
  //cout<< "alpha: "<< lp->getAlpha()<< endl;
  //cout<< "reversed? "<< lp->getReversed()<< endl;
  //cout<< "UpsideDown? "<< lp->getUpsideDown()<< endl;
  //cout<< "Quadrant" << lp->getQuadrant()<< endl;
  cout<< "cost " << lp->cost()<< endl;
}

void Problem:: checkQgsfeatureID(const bool& is_initial,unordered_map<int, int>& my_solution_prev){
  cout<< "check starts******************"<<endl;
  QSet<int> QgsFeatureIDS;
  int coutL = 0;
  int label;
  unordered_map<int, int> solution; 
  unordered_map<int, unordered_set<int>> candidates;
  //comparison with previous solution 
  if(!is_initial){
    if(numF_prev != nbft){
      cout<< "numF_prev: "<<numF_prev <<endl;  
      cout<< "nbft: "<<nbft <<endl;  
    };
    assert(numF_prev == nbft);
    for ( int i = 0; i < nbft; i++ ){
      coutL+=featNbLp[i];
      solution.insert({mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id(),sol->s[i]-featStartId[i]});
      candidates[mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id()];
      for ( int j = 0; j < featNbLp[i]; j++ ){  
        label = featStartId[i] + j;
        candidates[mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id()].insert(label);
        QgsFeatureIDS.insert(mLabelPositions.at(label)->getFeaturePart()->feature()->id());
      }
    }
    QSet<int>::const_iterator i;
    for (i = QgsFeatureIDS.begin(); i != QgsFeatureIDS.end(); ++i){
        if(!QgsFeatureIDS_prev.contains(*i)){
          cout<< "the current feature IDs:"<<endl;
          print(QgsFeatureIDS);
          cout<< "the previous feature IDs:"<<endl;
          print(QgsFeatureIDS_prev);
          cout<< "no "<< *i<< " in previous IDS"<< endl;
        }
        assert(QgsFeatureIDS_prev.contains(*i)); 
    }
    for (i = QgsFeatureIDS_prev.begin(); i != QgsFeatureIDS_prev.end(); ++i){
      if(!QgsFeatureIDS.contains(*i)){
          cout<< "the current feature IDs:"<<endl;
          print(QgsFeatureIDS);
          cout<< "the previous feature IDs:"<<endl;
          print(QgsFeatureIDS_prev);
          cout<< "no "<< *i<< " in current IDS"<< endl;

      }
      assert(QgsFeatureIDS.contains(*i)); 
    }
    if(coutL != numL_prev){
      cout<< "numL_prev: "<< numL_prev <<endl;  
      cout<< "current number: "<< coutL <<endl;
    }
    cout<< "CHECK THE SOLUTION"<<endl;
    int same = 0;
    int diff = 0;
    for(int i =0; i < nbft; i++){
      int qgsID = mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id();
      if(candidates_prev[qgsID].size() == candidates[qgsID].size()){
        if(my_solution_prev[qgsID] != solution[qgsID] && solution[qgsID] >= 0){
          printLabelPosition(mLabelPositions.at(featStartId[i]+solution[qgsID]));
          cout<<"the previous one is "<<  my_solution_prev[qgsID]<< endl;
        }
        else same++;
       /* if(solution[qgsID] <0){
          if(solution_prev[qgsID] >=0){
            printLabelPosition(mLabelPositions.at(featStartId[i]));
            QgsLabelFeature* qgsF =mLabelPositions.at(featStartId[i])->getFeaturePart()->feature();
            cout<<qgsF->labelText().toUtf8().constData() << endl;
          }
          assert(solution_prev[qgsID] < 0);
        }
        if(solution_prev[qgsID] <0){
          if(solution[qgsID] >=0){
            printLabelPosition(mLabelPositions.at(featStartId[i]+solution[qgsID]));
            QgsLabelFeature* qgsF =mLabelPositions.at(featStartId[i]+solution[qgsID])->getFeaturePart()->feature();
            cout<<qgsF->labelText().toUtf8().constData() << endl;
          }
          assert(solution[qgsID] < 0);
        }
        */
        if(candidates_prev[qgsID].size() == candidates[qgsID].size()){
          if(!(my_solution_prev[qgsID] == solution[qgsID])&&my_solution_prev[qgsID] >= 0 && solution[qgsID]>=0){
            diff++;
            cout<< "S: " << solution[qgsID]<< endl;
            cout<< "S_prev: " << my_solution_prev[qgsID]<< endl;
          }
          //assert((solution_prev[qgsID] == solution[qgsID])||(solution_prev[qgsID] < 0 || solution[qgsID]<0) );
        }
      }
    }
    cout<< "same "<< same << endl;
    cout<< "diff "<< diff<< endl;
    for(int i =0; i < nbft; i++){
      int qgsID = mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id();
      if(candidates_prev[qgsID].size() != candidates[qgsID].size()){
        QgsLabelFeature* qgsF = mLabelPositions.at(featStartId[i])->getFeaturePart()->feature();
        cout<<qgsF->labelText().toUtf8().constData() << endl;
      }
      //assert(candidates_prev[qgsID].size() == candidates[qgsID].size());
    }
  } 
  //assert(coutL == numL_prev);
  // store new solution
  QgsFeatureIDS_prev.clear();
  my_solution_prev.clear();
  candidates_prev.clear();
  //initial = false;
  numF_prev = nbft;
  coutL = 0;
  int fID_first = -1;
  int fID; 
    for (int  i = 0; i < nbft; i++ ){
      fID_first = -1;
      coutL+=featNbLp[i];
      my_solution_prev.insert({mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id(),sol->s[i]-featStartId[i]});
      candidates_prev[mLabelPositions.at(featStartId[i])->getFeaturePart()->feature()->id()];
      for ( int j = 0; j < featNbLp[i]; j++ ){ 
        label = featStartId[i] + j;
        fID =  mLabelPositions.at(label)->getFeaturePart()->feature()->id();
        candidates_prev[fID].insert(label);
        if(fID_first == -1) fID_first = fID;
        if(fID != fID_first){
          cout<< "fID: "<< fID<< endl;
          cout<< "fID_first: "<< fID_first<< endl;
        }
        assert(fID == fID_first);
        QgsFeatureIDS_prev.insert(fID);
      }
    }
  numL_prev = coutL;
  for(int i = 0; i < nbft; i++){
    cout<< featStartId[i]<< endl;    
    cout<< featNbLp[i]<< endl;
    for ( int j = 0; j < featNbLp[i]; j++ ){ 
              label = featStartId[i] + j;
              LabelPosition* LP = mLabelPositions.at(label);
              cout<< "getID: "<< LP->getId()<< endl;
              cout<< "label: " << label<< endl;

    }
  }
  cout<< "check ends******************"<<endl;
}
void Problem::checkLabelID(){
  int label;
  int label2;
  pal::LabelPosition::Quadrant quadrant;
  set<pal::LabelPosition::Quadrant> diff;
  for(int i = 0; i< nbft; i++){
    diff.clear();
    for(int j =0 ; j <featNbLp[i]; j++ ){ 
       label = featStartId[i] + j;
        quadrant = mLabelPositions.at(label)->getQuadrant();
        if(diff.find(quadrant)!= diff.end()){
              for(int j =0 ; j <featNbLp[i]; j++ ){ 
                cout<< "Quadrant is not unique!";
                label2 = featStartId[i] + j;
                printLabelPosition(mLabelPositions.at(label2));
              }
        }
        assert(diff.find(quadrant)== diff.end());
        diff.insert(quadrant);
    }
  }
}
void printFeature(QgsLabelFeature* qgsF){
  cout<< "Priority: "<< qgsF->priority()<< endl;
  cout<< "ID: "<< qgsF->id()<< endl;
  cout<<"Text:"<< qgsF->labelText().toUtf8().constData() << endl;

}
// priority for feature (used to define inativecost)and cost for label.
// same feature different label cost? (yes but very similiar check function createCandidate)
// there is no 4/8-position preferences for poit labeling
//TODO: define preference ourselves
// the cost is sorted by cost growing
void Problem::printCost(){
    int label;
     QgsLabelFeature* qgsF;
    for(int i = 0; i< nbft; i++){
      qgsF = mLabelPositions.at( featStartId[i])->getFeaturePart()->feature();
      //print feature info
      printFeature(qgsF);
      cout<< "has " << featNbLp[i] << "labels: "<< endl;
      for(int j =0 ; j <featNbLp[i]; j++ ){
        label = featStartId[i] + j;
        // print labelposition;
        printLabelPosition(mLabelPositions.at(label));
      }
    }
}
//-------------------Assertions for understanding QGIS-------------------------------
