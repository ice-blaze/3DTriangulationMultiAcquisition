#ifndef NODE_H
#define NODE_H

typedef int T;
const int K = 3;

class Node
{
    public:
        Node();
        Node(T _x, T _y, T _z);//TODO =0
        virtual ~Node();

        T key[K];

        Node* l;
        Node* r;

        static void add(Node** _root, Node* _element, int prof=0);
        static bool report(T _key[K], T _min[K], T _max[K]);
        static void search(Node* _root, T _min[K], T _max[K], int prof=0);

        T x();
        T y();
        T z();

        void x(T _x);
        void y(T _y);
        void z(T _z);
    protected:
    private:
};

#endif // NODE_H
