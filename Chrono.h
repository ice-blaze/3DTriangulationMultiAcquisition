#ifndef CHRONO_H
#define CHRONO_H

class Chrono
{
    public:
        Chrono();
        virtual ~Chrono();

        void restart();
        void printTime();
        int getTimeSeconde();

    protected:
    private:
        int start;
};

#endif // CHRONO_H
