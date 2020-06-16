/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <fstream>
#include <vector>
#include <iostream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "tmatrix.h"
#include "trackball.h"

using namespace std;

/*****Project 3 part*****/

#define	NPN 64
#define NMESH  100
#define DM  ((float) (1.0/(NMESH-1.0)))
#define NPIX  512
#define SCALE 4.0

/*****/

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;
float angle = 0.0f;
int max_num;
int single = 0;
double zoom = 1.0;
double translation[2] = {0, 0};

/*****The project 3 part*****/

int    iframe = 0;
int    Npat = 32;
int    alpha = (0.12 * 255);
float  sa;
float  tmax = NPIX / (SCALE * NPN);
float  dmax = SCALE / NPIX;

/*****/

struct jitter_struct {
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = { {0.0, 0.0} };
jitter_struct ji16[16] = { {0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125},
						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375},
						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625},
						  {0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly;

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);
void Animate();
/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
	char *progname;
	int num = 1;
	FILE *this_file;

	progname = argv[0];
	const char* xy = "../quadmesh_2D/output.ply";
	const char* mine = "../quadmesh_2D/sy/minecraft_sy.ply";

	this_file = fopen(mine, "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	mat_ident(rotmat);

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Geometric Modeling");
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutIdleFunc(Animate);
	glutMainLoop();
	poly->finalize();  // finalize everything

	return 0;    /* ANSI C requires main to return int. */
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0;
		col[1] = 1.0-percentage*3.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = 1.0;
		col[1] = percentage*3.0-1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 3.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file)
{
	int i, j;
	int elem_count;
	char *elem_name;

	/*** Read in the original PLY object ***/
	in_ply = read_ply(file);

	for (i = 0; i < in_ply->num_elem_types; i++) {

		/* prepare to read the i'th list of elements */
		elem_name = setup_element_read_ply(in_ply, i, &elem_count);

		if (equal_strings("vertex", elem_name)) {

			/* create a vertex list to hold all the vertices */
			nverts = max_verts = elem_count;
			vlist = new Vertex *[nverts];

			/* set up for getting vertex elements */

			setup_property_ply(in_ply, &vert_props[0]);
			setup_property_ply(in_ply, &vert_props[1]);
			setup_property_ply(in_ply, &vert_props[2]);
			setup_property_ply(in_ply, &vert_props[3]);
			setup_property_ply(in_ply, &vert_props[4]);
			setup_property_ply(in_ply, &vert_props[5]);
			setup_property_ply(in_ply, &vert_props[6]);

			vert_other = get_other_properties_ply(in_ply,
				offsetof(Vertex_io, other_props));

			/* grab all the vertex elements */
			for (j = 0; j < nverts; j++) {
				Vertex_io vert;
				get_element_ply(in_ply, (void *)&vert);

				/* copy info from the "vert" structure */
				vlist[j] = new Vertex(vert.x, vert.y, vert.z);
				vlist[j]->vx = vert.vx;
				vlist[j]->vy = vert.vy;
				vlist[j]->vz = vert.vz;

				vlist[j]->scalar = vert.s;

				vlist[j]->other_props = vert.other_props;
			}
		}
		else if (equal_strings("face", elem_name)) {

			/* create a list to hold all the face elements */
			nquads = max_quads = elem_count;
			qlist = new Quad *[nquads];

			/* set up for getting face elements */
			setup_property_ply(in_ply, &face_props[0]);
			face_other = get_other_properties_ply(in_ply, offsetof(Face_io, other_props));

			/* grab all the face elements */
			for (j = 0; j < elem_count; j++) {
				Face_io face;
				get_element_ply(in_ply, (void *)&face);

				if (face.nverts != 4) {
					fprintf(stderr, "Face has %d vertices (should be four).\n",
						face.nverts);
					exit(-1);
				}

				/* copy info from the "face" structure */
				qlist[j] = new Quad;
				qlist[j]->nverts = 4;
				qlist[j]->verts[0] = (Vertex *)face.verts[0];
				qlist[j]->verts[1] = (Vertex *)face.verts[1];
				qlist[j]->verts[2] = (Vertex *)face.verts[2];
				qlist[j]->verts[3] = (Vertex *)face.verts[3];
				qlist[j]->other_props = face.other_props;
			}
		}
		else
			get_other_element_ply(in_ply);
	}

	/* close the file */
	close_ply(in_ply);

	/* fix up vertex pointers in quads */
	for (i = 0; i < nquads; i++) {
		qlist[i]->verts[0] = vlist[(int)qlist[i]->verts[0]];
		qlist[i]->verts[1] = vlist[(int)qlist[i]->verts[1]];
		qlist[i]->verts[2] = vlist[(int)qlist[i]->verts[2]];
		qlist[i]->verts[3] = vlist[(int)qlist[i]->verts[3]];
	}

	/* get rid of quads that use the same vertex more than once */

	for (i = nquads - 1; i >= 0; i--) {

		Quad *quad = qlist[i];
		Vertex *v0 = quad->verts[0];
		Vertex *v1 = quad->verts[1];
		Vertex *v2 = quad->verts[2];
		Vertex *v3 = quad->verts[3];

		if (v0 == v1 || v1 == v2 || v2 == v3 || v3 == v0) {
			free(qlist[i]);
			nquads--;
			qlist[i] = qlist[nquads];
		}
	}
}


/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
	int i;
	PlyFile *ply;
	char **elist;
	int num_elem_types;

	/*** Write out the transformed PLY object ***/

	elist = get_element_list_ply(in_ply, &num_elem_types);
	ply = write_ply(file, num_elem_types, elist, in_ply->file_type);

	/* describe what properties go into the vertex elements */

	describe_element_ply(ply, "vertex", nverts);
	describe_property_ply(ply, &vert_props[0]);
	describe_property_ply(ply, &vert_props[1]);
	describe_property_ply(ply, &vert_props[2]);
	//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

	describe_element_ply(ply, "face", nquads);
	describe_property_ply(ply, &face_props[0]);

	//  describe_other_properties_ply (ply, face_other,
	//                                offsetof(Face_io,other_props));

	//  describe_other_elements_ply (ply, in_ply->other_elems);

	copy_comments_ply(ply, in_ply);
	char mm[1024];
	sprintf(mm, "modified by learnply");
	//  append_comment_ply (ply, "modified by simvizply %f");
	append_comment_ply(ply, mm);
	copy_obj_info_ply(ply, in_ply);

	header_complete_ply(ply);

	/* set up and write the vertex elements */
	put_element_setup_ply(ply, "vertex");
	for (i = 0; i < nverts; i++) {
		Vertex_io vert;

		/* copy info to the "vert" structure */
		vert.x = vlist[i]->x;
		vert.y = vlist[i]->y;
		vert.z = vlist[i]->z;
		vert.other_props = vlist[i]->other_props;

		put_element_ply(ply, (void *)&vert);
	}

	/* index all the vertices */
	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	/* set up and write the face elements */
	put_element_setup_ply(ply, "face");

	Face_io face;
	face.verts = new int[4];

	for (i = 0; i < nquads; i++) {

		/* copy info to the "face" structure */
		face.nverts = 4;
		face.verts[0] = qlist[i]->verts[0]->index;
		face.verts[1] = qlist[i]->verts[1]->index;
		face.verts[2] = qlist[i]->verts[2]->index;
		face.verts[3] = qlist[i]->verts[3]->index;
		face.other_props = qlist[i]->other_props;

		put_element_ply(ply, (void *)&face);
	}
	put_other_elements_ply(ply);

	close_ply(ply);
	free_ply(ply);
}

void Polyhedron::initialize() {
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize() {
	int i;

	for (i = 0; i < nquads; i++) {
		free(qlist[i]->other_props);
		free(qlist[i]);
	}
	for (i = 0; i < nedges; i++) {
		free(elist[i]->quads);
		free(elist[i]);
	}
	for (i = 0; i < nverts; i++) {
		free(vlist[i]->quads);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(qlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/

Quad *Polyhedron::find_common_edge(Quad *f1, Vertex *v1, Vertex *v2)
{
	int i, j;
	Quad *f2;
	Quad *adjacent = NULL;

	/* look through all faces of the first vertex */

	for (i = 0; i < v1->nquads; i++) {
		f2 = v1->quads[i];
		if (f2 == f1)
			continue;
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < f2->nverts; j++) {

			/* look for a match */
			if (f2->verts[j] == v2) {

#if 0
				/* watch out for triple edges */

				if (adjacent != NULL) {

					fprintf(stderr, "model has triple edges\n");

					fprintf(stderr, "face 1: ");
					for (k = 0; k < f1->nverts; k++)
						fprintf(stderr, "%d ", f1->iverts[k]);
					fprintf(stderr, "\nface 2: ");
					for (k = 0; k < f2->nverts; k++)
						fprintf(stderr, "%d ", f2->iverts[k]);
					fprintf(stderr, "\nface 3: ");
					for (k = 0; k < adjacent->nverts; k++)
						fprintf(stderr, "%d ", adjacent->iverts[k]);
					fprintf(stderr, "\n");

				}

				/* if we've got a match, remember this face */
				adjacent = f2;
#endif

#if 1
				/* if we've got a match, return this face */
				return (f2);
#endif

			}
		}
	}

	return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
	int i, j;
	Quad *f;

	/* make sure there is enough room for a new edge */

	if (nedges >= max_edges) {

		max_edges += 100;
		Edge **list = new Edge *[max_edges];

		/* copy the old list to the new one */
		for (i = 0; i < nedges; i++)
			list[i] = elist[i];

		/* replace list */
		free(elist);
		elist = list;
	}

	/* create the edge */

	elist[nedges] = new Edge;
	Edge *e = elist[nedges];
	e->index = nedges;
	e->verts[0] = v1;
	e->verts[1] = v2;
	nedges++;

	/* count all quads that will share the edge, and do this */
	/* by looking through all faces of the first vertex */

	e->nquads = 0;

	for (i = 0; i < v1->nquads; i++) {
		f = v1->quads[i];
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 4; j++) {
			/* look for a match */
			if (f->verts[j] == v2) {
				e->nquads++;
				break;
			}
		}
	}

	/* make room for the face pointers (at least two) */
	if (e->nquads < 2)
		e->quads = new Quad *[2];
	else
		e->quads = new Quad *[e->nquads];

	/* create pointers from edges to faces and vice-versa */

	e->nquads = 0; /* start this out at zero again for creating ptrs to quads */

	for (i = 0; i < v1->nquads; i++) {

		f = v1->quads[i];

		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 4; j++)
			if (f->verts[j] == v2) {

				e->quads[e->nquads] = f;
				e->nquads++;

				if (f->verts[(j + 1) % 4] == v1)
					f->edges[j] = e;
				else if (f->verts[(j + 2) % 4] == v1)
					f->edges[(j + 2) % 4] = e;
				else if (f->verts[(j + 3) % 4] == v1)
					f->edges[(j + 3) % 4] = e;
				else {
					fprintf(stderr, "Non-recoverable inconsistancy in create_edge()\n");
					exit(-1);
				}

				break;  /* we'll only find one instance of v2 */
			}

	}
}


/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
	int i, j;
	Quad *f;
	Vertex *v1, *v2;
	double count = 0;

	/* count up how many edges we may require */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++) {
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			Quad *result = find_common_edge(f, v1, v2);
			if (result)
				count += 0.5;
			else
				count += 1;
		}
	}

	/*
	printf ("counted %f edges\n", count);
	*/

	/* create space for edge list */

	max_edges = (int)(count + 10);  /* leave some room for expansion */
	elist = new Edge *[max_edges];
	nedges = 0;

	/* zero out all the pointers from faces to edges */

	for (i = 0; i < nquads; i++)
		for (j = 0; j < 4; j++)
			qlist[i]->edges[j] = NULL;

	/* create all the edges by examining all the quads */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < 4; j++) {
			/* skip over edges that we've already created */
			if (f->edges[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			create_edge(v1, v2);
		}
	}
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_quad_ptrs()
{
	int i, j;
	Quad *f;
	Vertex *v;

	/* zero the count of number of pointers to faces */

	for (i = 0; i < nverts; i++)
		vlist[i]->max_quads = 0;

	/* first just count all the face pointers needed for each vertex */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++)
			f->verts[j]->max_quads++;
	}

	/* allocate memory for face pointers of vertices */

	for (i = 0; i < nverts; i++) {
		vlist[i]->quads = (Quad **)
			malloc(sizeof(Quad *) * vlist[i]->max_quads);
		vlist[i]->nquads = 0;
	}

	/* now actually create the face pointers */

	for (i = 0; i < nquads; i++) {
		f = qlist[i];
		for (j = 0; j < f->nverts; j++) {
			v = f->verts[j];
			v->quads[v->nquads] = f;
			v->nquads++;
		}
	}
}


/******************************************************************************
Find the other quad that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Quad *Polyhedron::other_quad(Edge *edge, Quad *quad)
{
	/* search for any other quad */
	for (int i = 0; i < edge->nquads; i++)
		if (edge->quads[i] != quad)
			return (edge->quads[i]);

	/* there is no such other quad if we get here */
	return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_quad_ptrs(Vertex *v)
{
	int i, j;
	Quad *f;
	Quad *fnext;
	int nf;
	int vindex;
	int boundary;
	int count;

	nf = v->nquads;
	f = v->quads[0];

	/* go backwards (clockwise) around faces that surround a vertex */
	/* to find out if we reach a boundary */

	boundary = 0;

	for (i = 1; i <= nf; i++) {

		/* find reference to v in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[j] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #1\n");
			exit(-1);
		}

		/* corresponding face is the previous one around v */
		fnext = other_quad(f->edges[vindex], f);

		/* see if we've reached a boundary, and if so then place the */
		/* current face in the first position of the vertice's face list */

		if (fnext == NULL) {
			/* find reference to f in v */
			for (j = 0; j < v->nquads; j++)
				if (v->quads[j] == f) {
					v->quads[j] = v->quads[0];
					v->quads[0] = f;
					break;
				}
			boundary = 1;
			break;
		}

		f = fnext;
	}

	/* now walk around the faces in the forward direction and place */
	/* them in order */

	f = v->quads[0];
	count = 0;

	for (i = 1; i < nf; i++) {

		/* find reference to vertex in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[(j + 1) % f->nverts] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #2\n");
			exit(-1);
		}

		/* corresponding face is next one around v */
		fnext = other_quad(f->edges[vindex], f);

		/* break out of loop if we've reached a boundary */
		count = i;
		if (fnext == NULL) {
			break;
		}

		/* swap the next face into its proper place in the face list */
		for (j = 0; j < v->nquads; j++)
			if (v->quads[j] == fnext) {
				v->quads[j] = v->quads[i];
				v->quads[i] = fnext;
				break;
			}

		f = fnext;
	}
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Quad *f, Vertex *v)
{
	int j;
	int vindex = -1;

	for (j = 0; j < f->nverts; j++)
		if (f->verts[j] == v) {
			vindex = j;
			break;
		}

	return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_pointers()
{
	int i;

	/* index the vertices and quads */

	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	for (i = 0; i < nquads; i++)
		qlist[i]->index = i;

	/* create pointers from vertices to quads */
	vertex_to_quad_ptrs();

	/* make edges */
	create_edges();


	/* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
		order_vertex_to_quad_ptrs(vlist[i]);

	}
	/* index the edges */

	for (i = 0; i < nedges; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
		elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_sphere()
{
	unsigned int i;
	icVector3 min, max;

	for (i = 0; i < nverts; i++) {
		if (i == 0) {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else {
			if (vlist[i]->x < min.entry[0])
				min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
				max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
				min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
				max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
				min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
				max.entry[2] = vlist[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i = 0; i < nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1 - v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2, v3;
	Quad *temp_q;
	double edge_length[4];

	area = 0.0;
	for (i = 0; i < nquads; i++) {
		for (j = 0; j < 4; j++)
			edge_length[j] = qlist[i]->edges[j]->length;

		icVector3 d1, d2;
		d1.set(qlist[i]->verts[0]->x, qlist[i]->verts[0]->y, qlist[i]->verts[0]->z);
		d2.set(qlist[i]->verts[2]->x, qlist[i]->verts[2]->y, qlist[i]->verts[2]->z);
		double dia_length = length(d1 - d2);

		double temp_s1 = (edge_length[0] + edge_length[1] + dia_length) / 2.0;
		double temp_s2 = (edge_length[2] + edge_length[3] + dia_length) / 2.0;
		qlist[i]->area = sqrt(temp_s1*(temp_s1 - edge_length[0])*(temp_s1 - edge_length[1])*(temp_s1 - dia_length)) +
			sqrt(temp_s2*(temp_s2 - edge_length[2])*(temp_s2 - edge_length[3])*(temp_s2 - dia_length));

		area += qlist[i]->area;
		temp_q = qlist[i];
		v1.set(vlist[qlist[i]->verts[0]->index]->x, vlist[qlist[i]->verts[0]->index]->y, vlist[qlist[i]->verts[0]->index]->z);
		v2.set(vlist[qlist[i]->verts[1]->index]->x, vlist[qlist[i]->verts[1]->index]->y, vlist[qlist[i]->verts[1]->index]->z);
		v0.set(vlist[qlist[i]->verts[2]->index]->x, vlist[qlist[i]->verts[2]->index]->y, vlist[qlist[i]->verts[2]->index]->z);
		qlist[i]->normal = cross(v0 - v1, v2 - v1);
		normalize(qlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i = 0; i < nquads; i++) {
		icVector3 cent(vlist[qlist[i]->verts[0]->index]->x, vlist[qlist[i]->verts[0]->index]->y, vlist[qlist[i]->verts[0]->index]->z);
		signedvolume += dot(test - cent, qlist[i]->normal)*qlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume < 0)
		orientation = 0;
	else {
		orientation = 1;
		for (i = 0; i < nquads; i++)
			qlist[i]->normal *= -1.0;
	}
}

void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid) {
	unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

	if (sid >= eid)
		return;
	sort(A, B, C, sid, (sid + eid) / 2);
	sort(A, B, C, (sid + eid) / 2 + 1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	for (i = 0; i < eid - sid + 1; i++) {
		tempA[i] = A[i + sid];
		tempB[i] = B[i + sid];
		tempC[i] = C[i + sid];
	}
	current1 = sid;
	current2 = (sid + eid) / 2 + 1;
	current0 = sid;
	while ((current1 <= (sid + eid) / 2) && (current2 <= eid)) {
		if (tempA[current1 - sid] < tempA[current2 - sid]) {
			A[current0] = tempA[current1 - sid];
			B[current0] = tempB[current1 - sid];
			C[current0] = tempC[current1 - sid];
			current1++;
		}
		else if (tempA[current1 - sid] > tempA[current2 - sid]) {
			A[current0] = tempA[current2 - sid];
			B[current0] = tempB[current2 - sid];
			C[current0] = tempC[current2 - sid];
			current2++;
		}
		else {
			if (tempB[current1 - sid] < tempB[current2 - sid]) {
				A[current0] = tempA[current1 - sid];
				B[current0] = tempB[current1 - sid];
				C[current0] = tempC[current1 - sid];
				current1++;
			}
			else {
				A[current0] = tempA[current2 - sid];
				B[current0] = tempB[current2 - sid];
				C[current0] = tempC[current2 - sid];
				current2++;
			}
		}
		current0++;
	}
	if (current1 <= (sid + eid) / 2) {
		for (i = current1; i <= (sid + eid) / 2; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}
	if (current2 <= eid) {
		for (i = current2; i <= eid; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void init(void) {
	/* select clearing color */

	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 1;
		display();
		break;

	case '2':
		display_mode = 2;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	case 'x':
		switch (ACSIZE) {
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;

	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i = 0; i < 4; i++)
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i = 0; i < 4; i++)
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;

	case 'r':
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		display();
		break;

	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = nquads = 0;
	max_verts = max_quads = 50;

	vlist = new Vertex *[max_verts];
	qlist = new Quad *[max_quads];
}


void multmatrix(const Matrix m)
{
	int i, j, index = 0;

	GLfloat mat[16];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf(mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor * zoom, radius_factor*zoom, -radius_factor * zoom, radius_factor*zoom, 0.0, 40.0);
	else
		glFrustum(-radius_factor * zoom, radius_factor*zoom, -radius_factor, radius_factor, -1000, 1000);
		//gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(translation[0], translation[1], -3.0);
	multmatrix(rotmat);

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("Quad id = %d\n", seed_id);
	return seed_id;
}

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case -1:

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		glutPostRedisplay();
		break;

	case 3:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		glutPostRedisplay();
		break;

	case 2:

		zoom *= exp(2.0*(t - t_old));

		s_old = s;
		t_old = t;

		glutPostRedisplay();
		break;

	}
}

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		switch (mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_button = button;
				last_x = x;
				last_y = y;

				/*rotate*/
				if (button == GLUT_RIGHT_BUTTON)
				{
					mouse_mode = -1;
				}

				/*translate*/
				if (button == GLUT_LEFT_BUTTON)
				{
					mouse_mode = 3;
				}

				if (key == GLUT_ACTIVE_SHIFT)
				{
					mouse_mode = 2;
				}

			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
			GLint hits;
			GLint viewport[4];

			glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
			(void)glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			/*  create 5x5 pixel picking region near cursor location */
			gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			poly->seed = processHits(hits, selectBuf);
			glutPostRedisplay();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i = 0; i < poly->nquads; i++) {
		Quad *temp_q = poly->qlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
		for (j = 0; j < 4; j++) {
			Vertex *temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void Animate()
{
	angle += 0.05;
	//printf("called!\n");

	glutPostRedisplay();
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j, k;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	float colorX, colorY;
	double L = (poly->radius * 2) / 30;

	//Find max min to create the value sections
	double max = -99999;
	double min = 99999;
	double maxg = -99999;
	double ming = 99999;
	for (i = 0; i < this_poly->nverts; i++)
	{
		double thisv = sqrt((this_poly->vlist[i]->vx) * (this_poly->vlist[i]->vx)
			+ (this_poly->vlist[i]->vy) * (this_poly->vlist[i]->vy)
			+ (this_poly->vlist[i]->vz) * (this_poly->vlist[i]->vz));
		if (thisv < ming)
			ming = thisv;
		if (thisv > maxg)
			maxg = thisv;

		if (this_poly->vlist[i]->scalar < min)
			min = this_poly->vlist[i]->scalar;
		if (this_poly->vlist[i]->scalar > max)
			max = this_poly->vlist[i]->scalar;
	}

	//Define the number of sections
	//int numsec = sin(angle) * 100 + 100;
	int numsec = 40;

	//Create each section into the dynamic list
	double section = (max - min) / numsec;
	double *value = new double[numsec];
	for (i = 0; i < numsec; i++)
	{
		value[i] = min + section * i;
	}

	section = (maxg - ming) / numsec;
	double* valueg = new double[numsec];
	for (i = 0; i < numsec; i++)
	{
		valueg[i] = ming + section * i;
	}

	//Check the saddle point by checking each quads
	vector<Vertex> criticalList;
	vector<Vertex> criticalListg;
	vector<int> criticalSquad;
	vector<int> criticalSquadg;
	for (i = 0; i < this_poly->nquads; i++)
	{
		/*****The scalar part*****/
		Vertex* v1, * v2, * v3, * v4;
		v3 = this_poly->qlist[i]->verts[0];
		v4 = this_poly->qlist[i]->verts[1];
		v2 = this_poly->qlist[i]->verts[2];
		v1 = this_poly->qlist[i]->verts[3];

		double s1 = v1->scalar;
		double s2 = v2->scalar;
		double s3 = v3->scalar;
		double s4 = v4->scalar;

		double t = (s1 - s2) / (s1 - s2 - s3 + s4);
		double s = (s1 - s3) / (s1 - s2 - s3 + s4);

		//Check if the critical point inside the quad
		if (0 < t && t < 1 && 0 < s  && s < 1)
		{

			double ds = (v2->x - v1->x) * s;
			double dt = (v3->y - v1->y) * t;
			double x = v1->x + ds;
			double y = v1->y + dt;
			double z = 0;

			double sca = (1 - s) * (1 - t) * s1 + s * (1 - t) * s2 + (1 - s) * t * s3 + s * t * s4;

			Vertex newVert = Vertex(x, y, z);
			newVert.scalar = sca;
			//Store the scalar value
			criticalList.push_back(newVert);
			//Store the quad position in the list
			criticalSquad.push_back(i);
		}

		/*****The vxyz part*****/
		v3 = this_poly->qlist[i]->verts[0];
		v4 = this_poly->qlist[i]->verts[1];
		v2 = this_poly->qlist[i]->verts[2];
		v1 = this_poly->qlist[i]->verts[3];

		s1 = sqrt((v1->vx) * (v1->vx) 
			+ (v1->vy) * (v1->vy) 
			+ (v1->vz) * (v1->vz));
		s2 = sqrt((v2->vx) * (v2->vx)
			+ (v2->vy) * (v2->vy)
			+ (v2->vz) * (v2->vz));
		s3 = sqrt((v3->vx) * (v3->vx)
			+ (v3->vy) * (v3->vy)
			+ (v3->vz) * (v3->vz));
		s4 = sqrt((v4->vx) * (v4->vx)
			+ (v4->vy) * (v4->vy)
			+ (v4->vz) * (v4->vz));

		t = (s1 - s2) / (s1 - s2 - s3 + s4);
		s = (s1 - s3) / (s1 - s2 - s3 + s4);

		//Check if the critical point inside the quad
		if (0 < t && t < 1 && 0 < s && s < 1)
		{

			double ds = (v2->x - v1->x) * s;
			double dt = (v3->y - v1->y) * t;
			double x = v1->x + ds;
			double y = v1->y + dt;
			double z = 0;

			double sca = (1 - s) * (1 - t) * s1 + s * (1 - t) * s2 + (1 - s) * t * s3 + s * t * s4;

			Vertex newVert = Vertex(x, y, z);
			newVert.scalar = sca;
			//Store the scalar value
			criticalListg.push_back(newVert);
			//Store the quad position in the list
			criticalSquadg.push_back(i);
		}
	}

	//Check the vertex
	for (i = 0; i < this_poly->nverts; i++)
	{
		/*****The scalar part*****/
		Vertex* v = this_poly->vlist[i];
		Vertex* vu = NULL;
		Vertex* vr = NULL;
		Vertex* vd = NULL;
		Vertex* vl = NULL;
		int numNull = 4;
		int position = 0;
		if (this_poly->vlist[i]->nquads == 4)
		{
			vu = this_poly->vlist[i]->quads[0]->verts[1];
			vl = this_poly->vlist[i]->quads[0]->verts[3];
			vr = this_poly->vlist[i]->quads[1]->verts[2];
			vd = this_poly->vlist[i]->quads[2]->verts[3];
			double s = v->scalar;
			double su = vu->scalar;
			double sr = vr->scalar;
			double sd = vd->scalar;
			double sl = vl->scalar;
			if ((su > s && sd > s && sl < s && sr < s)
				|| (su < s && sd < s && sl > s && sr > s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalList.push_back(newVert);
				//Store the quad position in the list
				criticalSquad.push_back(-1);
			}
			else if ((su > s && sd > s && sl > s && sr > s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalList.push_back(newVert);
				//Store the quad position in the list
				criticalSquad.push_back(99);
			}
			else if ((su < s && sd < s && sl < s && sr < s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalList.push_back(newVert);
				//Store the quad position in the list
				criticalSquad.push_back(101);
			}
		}

		/*****The vxyz part*****/
		v = this_poly->vlist[i];
		vu = NULL;
		vr = NULL;
		vd = NULL;
		vl = NULL;
		numNull = 4;
		position = 0;
		if (this_poly->vlist[i]->nquads == 4)
		{
			vu = this_poly->vlist[i]->quads[0]->verts[1];
			vl = this_poly->vlist[i]->quads[0]->verts[3];
			vr = this_poly->vlist[i]->quads[1]->verts[2];
			vd = this_poly->vlist[i]->quads[2]->verts[3];
			double s = sqrt((v->vx) * (v->vx)
				+ (v->vy) * (v->vy)
				+ (v->vz) * (v->vz));
			double su = sqrt((vu->vx) * (vu->vx)
				+ (vu->vy) * (vu->vy)
				+ (vu->vz) * (vu->vz));
			double sr = sqrt((vr->vx) * (vr->vx)
				+ (vr->vy) * (vr->vy)
				+ (vr->vz) * (vr->vz));
			double sd = sqrt((vd->vx) * (vd->vx)
				+ (vd->vy) * (vd->vy)
				+ (vd->vz) * (vd->vz));
			double sl = sqrt((vl->vx) * (vl->vx)
				+ (vl->vy) * (vl->vy)
				+ (vl->vz) * (vl->vz));
			if ((su > s && sd > s && sl < s && sr < s)
				|| (su < s && sd < s && sl > s && sr > s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalListg.push_back(newVert);
				//Store the quad position in the list
				criticalSquadg.push_back(-1);
			}
			else if ((su > s && sd > s && sl > s && sr > s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalListg.push_back(newVert);
				//Store the quad position in the list
				criticalSquadg.push_back(99);
			}
			else if ((su < s && sd < s && sl < s && sr < s))
			{
				Vertex newVert = Vertex(v->x, v->y, v->z);
				newVert.scalar = s;
				//Store the scalar value
				criticalListg.push_back(newVert);
				//Store the quad position in the list
				criticalSquadg.push_back(101);
			}
		}
	}


	//printf("%d critical points found!\n", criticalList.size());

	bool ifDraw = true;
	for (i = 0; i < this_poly->nquads; i++) 
	{
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad *temp_q = this_poly->qlist[i];

		switch (display_mode) {
		case 0:
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_LIGHT1);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			}
			else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i == this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 1:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {
				Vertex *temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 2:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex *temp_v = temp_q->verts[j];
				
				colorX = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				colorY = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				glColor3f(colorX, colorY, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 3:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex *temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		case 4:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		case 5:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		case 6:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		case 7:
			glDisable(GL_LIGHTING);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glBegin(GL_POLYGON);
			for (j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);

			}
			glEnd();
			break;
		case 8:
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			for (j = 0; j < criticalList.size(); j++)
			{
				//Draw the contour point
				if (ifDraw)
				{
					//printf("Draw!\n");
					glPointSize(10.0);
					if(criticalSquad[j] == 99)
						glColor3f(1, 0, 1);
					else if(criticalSquad[j] == 101)
						glColor3f(0, 1, 0);
					else
						glColor3f(1, 0, 0);
					glBegin(GL_POINTS);
						glVertex3d(criticalList[j].x, criticalList[j].y, criticalList[j].z);
					glEnd();
				}

				//Find all the contour point and store in the list
				vector<Vertex> contourList;
				for (k = 1; k <= temp_q->nverts; k++)
				{
					Vertex* vert1 = temp_q->verts[k - 1];
					Vertex* vert2 = temp_q->verts[k];
					if (k == temp_q->nverts)
					{
						vert1 = temp_q->verts[k - 1];
						vert2 = temp_q->verts[0];
					}

					double v1 = vert1->scalar;
					double v2 = vert2->scalar;
					//Check two vertex
					if (((v1 - criticalList[j].scalar)) * ((v2 - criticalList[j].scalar)) <= 0)
					{
						//Find the relationship
						float ps = (criticalList[j].scalar - v1) / (v2 - v1);
						//Find coordinate of the point on contour
						double vert_x = vert1->x + (vert2->x - vert1->x) * ps;
						double vert_y = vert1->y + (vert2->y - vert1->y) * ps;
						double vert_z = vert1->z + (vert2->z - vert1->z) * ps;
						//Create Vertex and store in the list
						Vertex newContour = Vertex(vert_x, vert_y, vert_z);
						contourList.push_back(newContour);
					}
				}

				//Check if the squad contains the critical point
				if (i == criticalSquad[j])
					glBegin(GL_LINES);
				else
					glBegin(GL_LINE_STRIP);
						glLineWidth(1.5);
						glColor3f(0, 1, 1);
						for (k = 0; k < contourList.size(); k++)
						{
							if (i == criticalSquad[j])
							{
								glVertex3d(criticalList[j].x, criticalList[j].y, criticalList[j].z);
							}
							glVertex3d(contourList[k].x, contourList[k].y, contourList[k].z);
						}
					glEnd();
			}
			ifDraw = false;
			break;
			case 9:
			

				if (single == 0) {
					double a, b,c,d,e,f;
					for (j = 0; j < criticalList.size(); j++)
					{
						//Draw the contour point
						if (ifDraw)
						{
							if (criticalSquad[j] == 101) {
								printf("max point(x,y): %f,%f\n", criticalList[j].x, criticalList[j].y);
							}
							else if (criticalSquad[j] == 99) {
								printf("min point(x,y): %f,%f\n", criticalList[j].x, criticalList[j].y);

							}
							else {
								printf("sadle point(x,y): %f,%f\n", criticalList[j].x, criticalList[j].y);
							}
						}
						max_num = 0;
						for (int z = 0; z < criticalList.size(); z++)
						{
							for (k = z + 1; k < criticalList.size(); k++)
							{

								if (criticalSquad[z] == 101 && criticalSquad[k] == 101)
								{
										if (criticalList[z].y > criticalList[k].y) {
											a = criticalList[z].x;
											criticalList[z].x = criticalList[k].x;
											criticalList[k].x = a;

											b = criticalList[z].y;
											criticalList[z].y = criticalList[k].y;
											criticalList[k].y = b;
										}
									max_num++;
								}
								if (criticalSquad[z] == 99 && criticalSquad[k] == 99)
								{
									if (criticalList[z].y > criticalList[k].y) {
										c = criticalList[z].x;
										criticalList[z].x = criticalList[k].x;
										criticalList[k].x = c;

										d = criticalList[z].y;
										criticalList[z].y = criticalList[k].y;
										criticalList[k].y = d;
									}
								}
								if (criticalSquad[z] == -1 && criticalSquad[k] == -1)
								{
									if (criticalList[z].y > criticalList[k].y) {
										e = criticalList[z].x;
										criticalList[z].x = criticalList[k].x;
										criticalList[k].x = e;

										f = criticalList[z].y;
										criticalList[z].y = criticalList[k].y;
										criticalList[k].y = f;
									}
								}
							}
						}
					}
					printf("%d\n", max_num);
					printf("-----------------contour tree-------------\n");
					for (int j = 0; j < criticalList.size(); j++)
					{
						if (ifDraw)
						{
							if (criticalSquad[j] == 101) {
								int temp_num = j + max_num - 1;
								printf("max point(x,y): %f,%f\n", criticalList[temp_num].x, criticalList[temp_num].y);
								break;
							}
						}
					}
					printf("      |      \n");
					int start = 0;
					for (int j = 0; j < criticalList.size(); j++)
					{
						if (ifDraw)
						{
							if (criticalSquad[j] != 101 && criticalSquad[j] != 99) {
								if(start < max_num - 1){
									printf("sadle point(x,y): %f,%f", criticalList[j].x, criticalList[j].y);
									printf("      -      ");
									for (int k = 0; k < criticalList.size(); k++)
									{
											if (criticalSquad[k] == 101) {
												int temp_num = k + start;
												printf("max point(x,y): %f,%f\n", criticalList[temp_num].x, criticalList[temp_num].y);
												break;
											}
									}
									start++;
								}
								else{
									printf("sadle point(x,y): %f,%f\n", criticalList[j].x, criticalList[j].y);
								}
								printf("                    |      \n");
							}
						}
					}
				
					for (int j = 0; j < criticalList.size(); j++)
					{
						if (ifDraw)
						{
							if (criticalSquad[j] == 99) {
								printf("min point(x,y): %f,%f\n", criticalList[j].x, criticalList[j].y);
							}
						}
					}

					single++;
				}
				glDisable(GL_LIGHTING);
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				for (j = 0; j < criticalList.size(); j++)
				{
					//Draw the contour point
					if (ifDraw)
					{
						//printf("Draw!\n");
						glPointSize(10.0);
						if (criticalSquad[j] == 99) {
							glColor3f(1, 0, 1);
							
						}
						else if (criticalSquad[j] == 101)
							glColor3f(0, 1, 0);
						else
							glColor3f(1, 0, 0);
						glBegin(GL_POINTS);
						glVertex3d(criticalList[j].x, criticalList[j].y, criticalList[j].z);
						glEnd();
					}
					if (i == criticalSquad[j])
						glBegin(GL_LINES);
					else
						glBegin(GL_LINE_STRIP);
					glLineWidth(1.5);
					glColor3f(1, 0, 0);
					for (k = 0; k < criticalList.size(); k++)
					{
						if (criticalSquad[k] != 99 && criticalSquad[k] != 101) {
							if (i == criticalSquad[j])
							{

								glVertex3d(criticalList[k].x, criticalList[k].y, criticalList[k].z);
							}
							glVertex3d(criticalList[k].x, criticalList[k].y, criticalList[k].z);
						}

					}
					glEnd();

					glBegin(GL_LINE_STRIP);
					glLineWidth(1.5);
					glColor3f(1, 0, 1);
					int num1 = 0, num2 = 0, num3 = 0;
					for (k = 0; k < criticalList.size(); k++)
					{
						if (criticalSquad[k] == 99) {
							if (i == criticalSquad[j])
							{

								glVertex3d(criticalList[k].x, criticalList[k].y, criticalList[k].z);
								glColor3f(1.,1.,1.);
								glRasterPos2f(0., 0.);
								glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'a');
							}
							glVertex3d(criticalList[num1].x, criticalList[num2].y, criticalList[num3].z);
							num1++;
							num2++;
							num3++;
						}

					}
					glEnd();


					num1 = 5;
					num2 = 5;
					num3 = 5;
					for (k = 0; k < criticalList.size(); k++)
					{
						if (criticalSquad[k] == 101) {
							glBegin(GL_LINE_STRIP);
							glLineWidth(1.5);
							glColor3f(0, 1, 0);
							if (i == criticalSquad[j])
							{

								glVertex3d(criticalList[k].x, criticalList[k].y, criticalList[k].z);
							}
							glVertex3d(criticalList[num1].x, criticalList[num2].y, criticalList[num3].z);
							num1++;
							num2++;
							num3++;
							glEnd();
						}

					}
				}

					
				
				


				ifDraw = false;
				break;



		}
	}


}


void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor(0.0, 0.0, 0.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv(GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix();
		switch (ACSIZE) {
		case 1:
			glTranslatef(ji1[jitter].x*2.0 / viewport[2], ji1[jitter].y*2.0 / viewport[3], 0.0);
			break;

		case 16:
			glTranslatef(ji16[jitter].x*2.0 / viewport[2], ji16[jitter].y*2.0 / viewport[3], 0.0);
			break;

		default:
			glTranslatef(ji1[jitter].x*2.0 / viewport[2], ji1[jitter].y*2.0 / viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);
		display_shape(GL_RENDER, poly);
		glPopMatrix();
		glAccum(GL_ACCUM, 1.0 / ACSIZE);
	}
	glAccum(GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}

void Polyhedron::average_normals()
{
	int i, j;

	for (i = 0; i < nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j = 0; j < vlist[i]->nquads; j++)
			vlist[i]->normal += vlist[i]->quads[j]->normal;
		normalize(vlist[i]->normal);
	}
}


