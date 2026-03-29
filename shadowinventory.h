// ShadowInventory.h
#pragma once
#include "CoreMinimal.h"
#include "InventoryInterface.h"
#include "ShadowInventory.generated.h"

UCLASS(Blueprintable)
class UShadowInventory : public UObject, public IInventoryInterface
{
    GENERATED_BODY()

protected:
    UPROPERTY()
    TMap<FName, int64> Storage;

public:
    virtual bool AddItem_Implementation(TScriptInterface<IItem> Item, int32 Quantity) override;
    virtual bool RemoveItem_Implementation(FName ItemId, int32 Quantity) override;
    virtual int32 GetQuantity_Implementation(FName ItemId) const override;
};
