#include "Expense.h"

namespace Expenses
{
    Expense::Expense(std::string title, float_precision cost, tm* date)
    {
        this->expenseCost = cost;
        this->expenseDate = date;
    }

    float_precision Expense::getExpenseCost()
    {
        return this->expenseCost;
    }

    void Expense::setExpenseCost(float_precision costToSet)
    {
        this->expenseCost = costToSet;
    }

    void Expense::getExpenseDate(tm* timeptr)
    {
        *timeptr = *this->expenseDate;
    }
    void Expense::setExpenseDate(tm* timeptr)
    {
        this->expenseDate = timeptr;
    }
}