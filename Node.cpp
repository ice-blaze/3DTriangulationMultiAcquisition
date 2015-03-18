#include "Node.h"
#include <cstdio>

Node::Node()
{
    x(0);
    y(0);
    z(0);
}

Node::Node(T _x, T _y, T _z)
{
    x(_x);
    y(_y);
    z(_z);
}

Node::~Node()
{
}

void Node::add(Node** _root, Node* _element, int deep)
{
    if (*_root == nullptr)
        (*_root) = _element;
    else if ((*_root)->key[deep%K] < _element->key[deep%K])
        add(&(*_root)->r, _element, deep+1);
    else
        add(&(*_root)->l, _element, deep+1);
}

bool Node::report(T _key[K], T _min[K], T _max[K])
{
    int i, isIn = 1;
    for (i = 0; i < K && isIn; i++)
        isIn = isIn && _min[i] <= _key[i] && _max[i] >= _key[i];
    if (isIn) // Rapporter le point
    {
        for (i = 0; i < K; i++)
            printf("%f ", _key[i]);
        printf("\n");
    }
    return isIn;
}

void Node::search(Node* _root, T _min[K], T _max[K], int deep)
{
    if (_root) // arbre non vide
    {
        if (_root->key[deep%K] >= _min[deep%K])
            search(_root->l, _min, _max, deep+1);
        report(_root->key, _min, _max);
        if (_root->key[deep%K] <= _max[deep%K])
            search(_root->r, _min, _max, deep+1);
    }
}


T Node::x(){return key[0];}
T Node::y(){return key[1];}
T Node::z(){return key[2];}

void Node::x(T _x){key[0] = _x;}
void Node::y(T _y){key[1] = _y;}
void Node::z(T _z){key[2] = _z;}
