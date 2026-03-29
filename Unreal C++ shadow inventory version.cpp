// IItem.h
#pragma once
#include "UObject/Interface.h"
#include "IItem.generated.h"

UINTERFACE(Blueprintable)
class UItem : public UInterface
{
    GENERATED_BODY()
};

class IItem
{
    GENERATED_BODY()

public:
    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    FName GetItemId() const;

    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    FText GetItemName() const;

    UFUNCTION(BlueprintNativeEvent, BlueprintCallable)
    int32 GetMaxStack() const;
};
