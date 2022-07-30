#include "Account.h"
using namespace Accounts;

Account::Account()
{
}

float_precision Account::getBalance()
{
    return this->balance;
}

void Account::setBalance(float_precision toSet)
{
    this->balance = toSet;
}
