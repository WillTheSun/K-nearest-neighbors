 //William Sun
//12-3-2014 
//Implementation of a Quad Tree to find the nearest neighbors
//of all points in the array

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define FIND_DISTANCE(a,b) sqrt((b->x - a->x)*(b->x - a->x)+ (b->y - a->y)*(b->y - a->y))

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct point {
	double x;
	double y;
	int index;
	int data;
} point;

typedef struct node {
	int numPoints;
	double length;
	point *center;

	point *nwPoint;
	point *nePoint;
	point *swPoint;
	point *sePoint;
	point **pointList;
	struct node *parent;
	struct node *nw;
	struct node *ne;
	struct node *sw;
	struct node *se;
} node;

point **createPointList(int n){
	return malloc(n * sizeof(point*));
}

double *createDoubleList(int k){
	double *List;
	List = malloc(k * sizeof(double));
		int l;
		for(l = 0; l < k; l++)
			List[l] = 100000000000000000;
	return List;
}

//SET POINT VALUES
point *setPoint(double x, double y, int index) {
	point *p = malloc(sizeof(point));
	p->x = x;
	p->y = y;
	p->index = index;
	return p;
}

// Create NODE
node *createNode(node *parent, int number, point *sw, point **points, double length) {
	node *n = malloc(sizeof(node));
	n->parent = parent;
	n->numPoints = number;
	n->length = length;
	n->nwPoint = setPoint((*sw).x, (*sw).y + length, -1);
	n->nePoint = setPoint((*sw).x + length, (*sw).y + length, -1);
	n->swPoint = setPoint((*sw).x, (*sw).y, -1);
	n->sePoint = setPoint((*sw).x + length, (*sw).y, -1);
	n->pointList = points;
	return n;
}

// does the circle intersect with a N?
int isInCircle(point *a, double radius, node *N) {
	double nx = (a->x + radius);
	double nx2 = (a->x - radius);
	double ny = a->y + radius;
	double ny2 = a->y - radius;
	if(FIND_DISTANCE(a, N->nwPoint) <= radius ||
		FIND_DISTANCE(a, N->nePoint) <= radius ||
		FIND_DISTANCE(a, N->swPoint) <= radius ||
		FIND_DISTANCE(a, N->sePoint) <= radius)
		return 1;
	return 0;
}

// Insert the point into Iz
void insertPoint(node *N, double *list, int k, point *a, int *iz) {
	int i,j;
	for(i=0; i<N->numPoints; i++){
		point *P = N->pointList[i];
		int index = P->index;
		double minX = 100000000000000;
		double minY = 100000000000000;
		double maxX = -100000000000000;
		double maxY = -100000000000000;
		double D = FIND_DISTANCE(a, N->pointList[i]);
		j=0;
		if(D == 0)
			continue;
		if (D >= list[j])
			break;
		else{
			while((D < list[j]) && (j < k-1)){
				if(D > list[j+1])
					list[j] = D;
				else{
					list[j] = list[j+1];
					iz[(a->index)*k + j] = iz[(a->index)*k +j+1];
					j++;
				}
			}
			iz[(a->index)*k + j] = N->pointList[i]->index;
			list[j] = D;
			minX = minY = D;
		}
	}
}

// Grabs every point frmo every node
void getPossiblePoints(node *N, double *List, int k, point *a, double radius, int *iz) {
	// Still a branch?
	if(N->nw != NULL){
		if(isInCircle(a, radius, N->nw))
		getPossiblePoints(N->nw, List, k, a, radius, iz);
	if(isInCircle(a, radius, N->ne))
		getPossiblePoints(N->ne, List, k, a, radius, iz);
	if(isInCircle(a, radius, N->sw))
		getPossiblePoints(N->sw, List, k, a, radius, iz);
	if(isInCircle(a, radius, N->se))
		getPossiblePoints(N->se, List, k, a, radius, iz);
	}
	else 
		insertPoint(N, List, k, a, iz);
}

// We have the largest N. Set up the next recursive function to 
// check all points within the sub leaves
void findNeighbors(node *N, point *a, int k, int *iz) {
	int i, j;
	double nwRadius = FIND_DISTANCE(a, N->nwPoint);
	double neRadius = FIND_DISTANCE(a, N->nePoint);
	double swRadius = FIND_DISTANCE(a, N->swPoint);
	double seRadius = FIND_DISTANCE(a, N->sePoint);
	double *List;
	List = createDoubleList(k);
	// First find the largest N
	while(isInCircle(a, neRadius, N) && isInCircle(a, nwRadius, N) &&
		isInCircle(a, swRadius, N) && isInCircle(a, seRadius, N) &&
		N->parent != NULL){
		N = N->parent;
		N->center = 0;
	}
	double radius = MAX(MAX(nwRadius, neRadius), MAX(swRadius, seRadius));
	getPossiblePoints(N, List, k, a, radius, iz);
	free(List);
}

//Recurses to Leaf level
void findLeaves(node *N, int k, int *iz) {
	int i;
	//Recursively search other trees
	if (N->nw != NULL){
		findLeaves(N->nw, k, iz);
		findLeaves(N->ne, k, iz);
		findLeaves(N->sw, k, iz);
		findLeaves(N->se, k, iz);
	}
	for(i=0; i< N->numPoints; i++) 
		findNeighbors(N, N->pointList[i], k, iz);	
}

void createTree(node *N, int k) {
	int NN = N->numPoints;
	int size[4] = {0,0,0,0};

	if(NN <= k)
		return;
	int i;
	double newL = (N->length)/2.0;
	point *sw = N->swPoint;
	point **swPointList = createPointList(NN);
	point **sePointList = createPointList(NN);
	point **nwPointList = createPointList(NN);
	point **nePointList = createPointList(NN);
	for(i = 0; i < NN; i++) {
		point *P = N->pointList[i];
		if(P->x < (sw->x + newL) && P->y < (sw->y + newL)){
			swPointList[size[0]] = P;
			size[0]++;
		}
		else if(P->x >= (sw->x+newL) && P->y < (sw->y + newL)){
			sePointList[size[1]] = P;
			size[1]++;
		}
		else if(P->x >= (sw->x + newL) && P->y >= (sw->y + newL)){
			nePointList[size[2]] = P;
			size[2]++;
		}
		else{
			nwPointList[size[3]] = P;
			size[3]++;
		}
	}
	point *p1 = setPoint(sw->x + newL, sw->y, -1);
	point *p2 = setPoint(sw->x + newL, sw->y + newL, -1);
	point *p3 = setPoint(sw->x, sw->y + newL, -1);

	N->sw = createNode(N, size[0], sw, swPointList, newL);
	N->se = createNode(N, size[1], p1, sePointList, newL);
	N->ne = createNode(N, size[2], p2, nePointList, newL);
	N->nw = createNode(N, size[3], p3, nwPointList, newL);
	createTree(N->sw, k);
	createTree(N->se, k);
	createTree(N->ne, k);
	createTree(N->nw, k);
}

void seek(double *a, int n, int k, int *iz) {
	int i;
	point **pointList = malloc(n * sizeof(point*));
	point *sw = (point*) malloc(sizeof(point));
	point *nw = (point*) malloc(sizeof(point));
	for(i = 0; i < n; i++)
		pointList[i] = setPoint(a[2 * i], a[2 * i + 1], i);
	double minX = pointList[0]->x, maxX = pointList[0]->x;
	double minY = pointList[0]->y, maxY = pointList[0]->y;
	for(i = 0; i < n; i++) {
		minX = MIN(pointList[i]->x, minX);
		maxX = MAX(pointList[i]->x, maxX);
		minY = MIN(pointList[i]->y, minY);
		maxY = MAX(pointList[i]->y, maxY);
	}
	sw = setPoint(minX, minY, -1);
	double length = MAX(maxX - minX, maxY - minY);
	node *root;
	root = (node*) malloc(sizeof(node));
	root = createNode(NULL, n, sw, pointList, length);
	createTree(root, k);
	findLeaves(root, k, iz);
	free(pointList);
	free(sw);
	free(root);
}

void seek_naive(double *a, int n, int k, int *iz){
	int i, j, l;
	point **pointList;
	pointList = malloc(n * sizeof(point*));
	for(i=0; i<n; i++)
		pointList[i] = setPoint(a[2*i], a[2*i+1], i);
	for(i = 0; i<n; i++) {
		double *List = createDoubleList(k);
		for(j=0; j<n; j++){
			l=0;
			if(i == j)
				continue;
			double D = FIND_DISTANCE(pointList[i], pointList[j]);
			if(D < List[l]) {
				while(D < List[l] && l < k-1) {
					if(D > List[l])
						List[l] = D;
					else {
						List[l] = List[l+1];
						iz[(pointList[i]->index)*k + l] = iz[(pointList[i]->index)*k + l + 1];
						l++;
					}
				}
				List[l] = D;
				iz[(pointList[i]->index)*k + l] = pointList[j]->index;
			}
		}
		free(List);
	}
}

int main() {
	int i,j;
	int k = 3;
	int n = 50;
	double *a = malloc(2*n*sizeof(double));
	int *iz = malloc(k*n*sizeof(int));
	for (i=0;i<2*n;i++){
		a[i] = n*i - i*i;
	}
	for(i=0;i<2*n;i=i+2){
		a[i] = -sqrt(i*i*i)/100;
	}
	seek(a, n, k, iz);
	for (i=0;i<n;i++){
		printf("Point (%lf,%lf):\n", a[i],a[2*i]);
		for (j=0;j<k;j++){
			int index = iz[i*k + j];
			printf("\tNeighbor %d: (%lf,%lf)\n", j, a[index],a[2*index]);
		}
	}
}
