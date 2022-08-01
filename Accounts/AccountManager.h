//#include <vector>
#ifndef ACCOUNTMANAGER_H_
#define ACCOUNTMANAGER_H_
#include "Account.h"
#include "Expenses.h"
using namespace Accounts;
using namespace Expenses;

namespace AccountManagers
{
    class AccountManager
    {
    public:
        void addAccount();
        Account getAccount(int index);

    private:
        std::vector<Account> AccountList;
    }
}

#endif // ACCOUNTMANAGER_H_